# Agent Guide

## Scope
- Repo combines a legacy Flask app (`app/`, Python 2.7) with the REST API
  (`app_v1/`, Python 3.6) and the tool runners in `v2/`, deployed to `/opt/v2`.
- Several Python versions coexist: 2.7 (legacy stack, `mirna_processing.py`),
  3.6 (`app_v1`, DMISO, most tool runners), 3.11 (miRmap 2 only).

## Python Version Guidance
- `v2/mirna_processing.py` runs under `python2.7` and MUST stay 2.7-compatible
  (no f-strings, no `exist_ok=True`, no `subprocess.run`).
- `v2/mirna_predicting.py` is invoked as `python3.6` (every tool except miRmap)
  and `python3.11` (miRmap only), so it must import cleanly under both.
- All of `v2/*.py` must be ASCII, or carry a PEP 263 encoding declaration on
  line 1–2 — Python 2.7 raises `SyntaxError` otherwise. Enforced by
  `v2/tests/test_source_encoding.py` and the `.githooks` pre-commit hook
  (`git config core.hooksPath .githooks`).
- miRmap is **miRmap 2**, pip-installed into python3.11 site-packages. miRmap 1.x
  is gone from the image. `mirna_predicting.py` strips `/opt/miRmap/src` from
  `sys.path` under Python 3 so the stale 1.x source cannot shadow it.
- DMISO runs via the `/usr/local/bin/dmiso` wrapper, which execs
  `python3.6 /opt/DMISO/DMISO-main/dmiso.py`.

## Containers & Docker
- Base image file: `isotar-base.Dockerfile` (slow apt installs, Python 3.6 and
  3.11 built from source, miRmap 2). Its tag must match the `FROM` in `Dockerfile`.
- Final image file: `Dockerfile` (tool setup, app setup).
- Build sequence:
  - `docker build -t frankfeng78/isotar-v2-base:0.2.x -f isotar-base.Dockerfile .`
  - `docker build -t frankfeng78/isotar-v2:"$(cat VERSION)" --build-arg VERSION="$(cat VERSION)" -f Dockerfile .`
- Releases are cut with `scripts/release.sh <major|minor|patch|X.Y.Z>`; it
  refuses to run with uncommitted *tracked* changes (untracked files are fine).

## Paths & Data
- 3' UTR datasets: `/opt/reference_files/<code>_<assembly>_3UTRs.fasta`
  (e.g. `hsa_HG19_only_3UTRs.fasta`, `hsa_HG38_only_3UTRs.fasta`), shipped from
  `v2/opt/reference_files/`.
- lncRNA datasets: `/opt/reference_files/<code>_<assembly>_lncRNAs.fasta`. These
  are NOT committed — build them with `scripts/build_lncrna_references.sh`.
- miRNA metadata: `/opt/resources/` (`ISOTAR_RESOURCES_DIR`). Human is the
  un-prefixed `mature_pre_mirna_ext.json`; other species are
  `<code>_mature_pre_mirna_ext.json`, resolved at runtime from the miRBase
  prefix in the miRNA id.
- TargetScan datasets live under `/opt/TargetScan/Datasets`, copied in by
  `ADD tools /opt/` from `tools/TargetScan/Datasets/` (gitignored, so not in the
  repo). TargetScan reads THESE, never `/opt/reference_files` — its 3' UTR
  boundaries differ, so swapping references changes what it predicts rather than
  extending it (on hsa-miR-21-5p only 33% of shared transcripts have the same
  UTR length; the two references agree on 70% of genes).
  - `Datasets/3utr/` is human and serves both `hg19` and `hg38`.
  - `Datasets/<genome>/3utr/` is per species, built by
    `scripts/build_targetscan_species_datasets.sh` (~463 MB download, ~4.3 GB on
    disk). Each is 64 parts split on gene boundaries — the shipped human
    splitter's fixed 37,212 lines only works because human is exactly 84 rows
    per gene (mouse is 52, fly 27, zebrafish 1).
  - `Datasets/bln_bins/` is human-only. `targetscan_70_BL_bins.pl` cannot be run
    in the image (needs `Statistics::Lite`, not installed), so other genomes get
    a flat bin table derived from script 1's output — see `_write_flat_bins`.
- TargetScan runs for `hg19`, `hg38`, `mmu`, `dme`, `dre` only, and the API
  returns 400 for the rest. It needs an identifier it can map back to RefSeq:
  worm keys every file on an internal numeric (`171590.0`) with no RefSeq or
  Ensembl equivalent, and `rno`/`cfa`/`mml`/`ptr`/`mdo` have no TargetScan
  release at all. The map is the `enst_refseq` table in
  `app_v1/reference_mapping.db`, keyed by assembly (`GRCh37`, `GRCh38`,
  `GRCm38`, `Release6`, `GRCz11`) — ENST, ENSMUST, FBtr and ENSDARG accessions
  from Ensembl BioMart, curated `refseq_mrna` only.
- Jobs live under `ISOTAR_JOB_DIR` (`/opt/out/jobs`), one directory per job.

## Job Execution
- A job is a directory: `job.json` (meta + status), `mirna.fa`, optional
  `targets.txt`, and `output/` (raw tool files + `progress.json`).
- `v2/mirna_processing.py` generates the miRNA FASTA, including isomiR variants
  (`-m <pos>:<from>|<to>` modifications, `-s left|right` precursor shift).
- `v2/mirna_predicting.py` runs tools based on `-t`, against the pool named by
  `-tt` (`gene` | `lncrna`), optionally narrowed by `-tf targets.txt`.
- Parallel branch uses multiprocessing `pool.map` with wrapper functions.
- `app_v1` runs the tools as two subprocess groups per pool (python3.6 tools,
  then python3.11 miRmap). One group failing must not discard the other's
  results — failures land in `job.json`'s `tool_errors`, and the job only fails
  when every group failed.
- Three workflows: `mir-target` (gene pool), `mir-lncrna` (lncRNA pool,
  TargetScan rejected), `mir-network` (both pools into `output/gene/` and
  `output/lncrna/`, then joined into a tripartite graph).
- Results are parsed lazily: the first `/result` call builds and caches
  `result.db` in the output directory.

## Common Issues
- DMISO wrapper must preserve args; call the script directly if needed.
- `subprocess.run` isn't available in Python 2; ensure a shim exists where used
  on the 2.7 path.
- PITA needs a working `RNAddG4` accessibility binary — if it produces no
  output, `dGopen` collapses and predictions over-predict by 20–40x. The parser
  raises rather than returning silently wrong results
  (`scripts/rebuild_RNAddG4.sh`). Each PITA run also needs a private cwd.
- TargetScan ignores the supplied FASTA and reads its own precomputed datasets,
  so it cannot run against the lncRNA pool, and it only supports the genomes it
  ships data for (see Paths & Data).
- TargetScan fails SILENTLY on a species mismatch. `targetscan_70.pl` skips every
  UTR row whose species id is absent from the miRNA file, so a wrong id gives an
  empty result and exit 0 — not an error. Two places must agree with the job's
  genome: the species list written by `targetscan_prep` (the shipped
  `miR_Family_Info` is TargetScan's VERTEBRATE set, so a fly seed can match it
  and come back as 9606 — the genome's own taxon is always unioned in), and the
  `species_id` row filter in `parseTargetScanResults`. Filtering mouse output at
  9606 returns 0 targets rather than erroring, because a mouse alignment
  contains human rows.
- Filenames carry the miRNA variant tag (`<mirna>,<variant>_<tool>_results.txt`);
  sanitize before handing them to tools that choke on the punctuation.

## Testing
- `v2/tests/` is plain `unittest`, no `__init__.py` (implicit namespace package).
  Run from the repo root by naming modules; `unittest discover` will not work:
  - `python3 -m unittest v2.tests.test_parse_result`
  - `python3 -m unittest $(ls v2/tests/test_*.py | sed 's|/|.|g; s|\.py$||')`
- Modules importing `app_v1` need Flask/Celery present, otherwise they error at
  import time.

## Expectations
- Prefer minimal changes, keep legacy behavior.
- Avoid upgrading libraries unless required for compatibility.
- Keep `CLAUDE.md` and `AGENTS.md` identical — they are two copies of this file.
