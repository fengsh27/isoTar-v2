# isoTar-v2

A microRNA (miRNA) target prediction platform that runs six computational tools in parallel and aggregates results into a searchable database. Supports canonical miRNA sequences as well as isomiR variants (modified, shifted, or truncated).

## Overview

isoTar-v2 accepts a miRNA identifier (or a raw custom sequence) plus optional
sequence modifications, generates the corresponding FASTA, and dispatches
prediction jobs to any combination of six tools. Jobs are processed
asynchronously via Celery. Results are stored in SQLite and exposed through a
REST API with pagination, sorting, and keyword filtering.

Three workflows share the same job machinery — see [Workflows](#workflows).

## Architecture

```
Client
  └── HTTP :8080
        └── Nginx (reverse proxy)
              ├── /api/v1/*  →  Gunicorn :5001  (app_v1 — Python 3.6)
              └── /*         →  Gunicorn :5000  (legacy app — Python 2.7)

app_v1 (Flask + Celery)
  └── run_job task
        ├── python2.7  mirna_processing.py   (FASTA generation)
        └── per target pool, two independent subprocess groups:
              ├── python3.6   mirna_predicting.py  (all tools except miRmap)
              └── python3.11  mirna_predicting.py  (miRmap 2 only)

Supervisor manages: nginx, rabbitmq, gunicorn (legacy), app_v1 gunicorn, celery worker
```

A job is a directory, not a database row: `$ISOTAR_JOB_DIR/<job_id>/` holds
`job.json` (metadata + status), `mirna.fa`, an optional `targets.txt`, and
`output/` with the raw per-tool files and `progress.json`. The two tool groups
run independently — if one fails the other's results are kept, and the failure
is recorded in `job.json` under `tool_errors`. A job is only marked `failed`
when *every* group failed.

Nothing is parsed until the first `/result` call, which builds and caches
`output/result.db` from the raw tool output.

## Prediction Tools

| Tool | Version | Method |
|------|---------|--------|
| miRanda | 3.3 | Seed-based thermodynamic scoring |
| miRmap | 2.x | Comprehensive thermodynamic model (pip-installed, runs under Python 3.11) |
| RNAhybrid | 2.1.2 | RNA–RNA interaction energy |
| PITA | v6 | Position-and-context-dependent ddG scoring |
| TargetScan | 7.0 | Seed match + conservation scoring |
| DMISO | — | Deep learning (TensorFlow/Keras) |

## Workflows

Every job carries a `workflow` field selecting which target sequence pool the
tools are scored against. It defaults to `mir-target`.

| `workflow` | Target pool | Notes |
|------------|-------------|-------|
| `mir-target` | gene 3′ UTRs | The default. One miRNA per job. |
| `mir-lncrna` | lncRNA transcripts | One miRNA per job. TargetScan is rejected with HTTP 400 — it ignores the supplied FASTA and reads its own precomputed 3′ UTR datasets. |
| `mir-network` | **both** pools | A list of miRNAs (max `ISOTAR_MAX_NETWORK_MIRNAS`) plus optional ceRNA `(gene, lncRNA)` pairs with correlation magnitudes, joined into a tripartite gene ↔ miRNA ↔ lncRNA graph and scored. |

A `mir-network` job runs each pool into its own output subdirectory
(`output/gene/`, `output/lncrna/`), each with its own `progress.json`; the API
merges them for `GET /jobs/<id>`. Because there is no single result set at the
top level, `GET /jobs/<id>/result` **requires** a `pool` parameter for these
jobs, and the graph itself is served by `GET /jobs/<id>/network`.

## Supported Genomes / Species

| Code | Species | TargetScan |
|------|---------|:----------:|
| `hg19` | Human (GRCh37) | yes |
| `hg38` | Human (GRCh38) | yes |
| `mmu` | House mouse (GRCm38) | yes |
| `rno` | Norway rat (RGSC6/rn6) | — |
| `dre` | Zebrafish (GRCz11) | yes |
| `dme` | Fruit fly (Release 6) | yes |
| `cel` | Roundworm (WBcel235) | yes |
| `cfa` | Dog (CanFam3.1) | — |
| `mdo` | Gray short-tailed opossum (MonDom5) | — |
| `mml` | Rhesus macaque (Mmul_8.0.1) | — |
| `ptr` | Chimpanzee (Pan_tro3.0) | — |

Five of the six tools run on every genome above. **TargetScan is the exception**:
it ignores the supplied FASTA and scores against its own precomputed 3′ UTR
alignments, so it only runs where TargetScan publishes a dataset whose
identifiers can be mapped back to RefSeq. Requesting it for any other genome
returns HTTP 400; rat, dog, opossum, macaque and chimpanzee have no TargetScan
release at all. Roundworm's TargetScan identifiers name a gene rather than a
transcript, so each worm hit is reported against every RefSeq transcript of its
gene (1.24 on average), as zebrafish's already are. Its per-species datasets are
not committed — build them with `scripts/build_targetscan_species_datasets.sh`
(~469 MB download, ~4.3 GB in the image).

3′ UTR references ship in the repo under `v2/opt/reference_files/` and land at
`/opt/reference_files/<code>_<assembly>_3UTRs.fasta` in the image. The lncRNA
references (`*_lncRNAs.fasta`, needed by the `mir-lncrna` and `mir-network`
workflows) are **not** committed — build them from the Ensembl ncrna dumps with
`scripts/build_lncrna_references.sh`.

## Versioning

The current version lives in the [`VERSION`](VERSION) file (Semantic Versioning)
and drives the `frankfeng78/isotar-v2` image tag, the `vX.Y.Z` git tag, and the
`/api/v1/version` endpoint. Changes are recorded in [`CHANGELOG.md`](CHANGELOG.md).

Cut a release (bumps `VERSION`, updates the changelog, commits and tags):

```bash
scripts/release.sh patch     # or: minor | major | 0.3.0
```

Check what a running container reports:

```bash
curl -s http://localhost:8080/api/v1/version
# {"name":"isotar-v2-backend","version":"0.2.27"}
```

## Docker Build

```bash
# Step 1 — base image (Ubuntu 16.04; Python 3.6 and 3.11 compiled from source,
# build tools, miRmap 2). Must match the FROM tag in Dockerfile.
# Only needed once or when base dependencies change
docker build -t frankfeng78/isotar-v2-base:0.2.3 -f isotar-base.Dockerfile .

# Step 2 — final production image (stamp the version into the image)
docker build -t frankfeng78/isotar-v2:"$(cat VERSION)" \
  --build-arg VERSION="$(cat VERSION)" -f Dockerfile .
```

## Running

```bash
cp .env.template .env        # fill in NCBI and gateway credentials (optional)
docker run -d \
  --name isotar \
  -p 8080:8080 \
  --env-file .env \
  -v /path/to/jobs:/opt/out/jobs \
  frankfeng78/isotar-v2:"$(cat VERSION)"
```

The API is available at `http://localhost:8080/api/v1/`.

`.env.template` lists every environment variable the server reads, with the
defaults and a note on what each one changes. The two optional feature blocks
are the NCBI E-utilities credentials and the Anthropic-compatible gateway
(LiteLLM) used by the literature verification; see
[Environment variables](#environment-variables). `.env` is git-ignored.

## REST API

All endpoints are under `/api/v1`. Bodies and responses are JSON unless noted.

### Validate targets

```
POST /api/v1/targets/validate
```

Read-only pre-flight check — creates no job state. Use it to tell a user which
of their targets are unknown before they submit a job (submission rejects the
whole job if any target is unresolved).

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `targets` | string[] | yes | Up to `ISOTAR_MAX_VALIDATE_TARGETS` identifiers |
| `genome` | string | no | Species code (default `hg38`). Requesting `Targetscan` for a genome it has no dataset for returns HTTP 400 — see [Supported Genomes](#supported-genomes--species) |
| `target_type` | string | no | `gene` (default) or `lncrna` |

`gene` checks HGNC symbols and RefSeq accessions against `gene_mapping`;
`lncrna` checks Ensembl / FlyBase / WormBase lncRNA transcript and gene IDs.

```json
{
  "genome": "hg38",
  "species": "hsa",
  "target_type": "gene",
  "results": [{ "target": "TP53", "valid": true, "matched_by": "symbol" }],
  "valid_count": 1,
  "invalid": []
}
```

`matched_by` is `symbol` or `accession` for `gene` targets, `transcript` or
`gene` for `lncrna` targets, and `null` when the target did not resolve.

---

### Submit a job

```
POST /api/v1/jobs
```

Common fields:

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `tools` | string[] | yes | Tools to run: `miRanda`, `miRmap`, `RNAhybrid`, `PITA`, `Targetscan`, `DMISO` |
| `workflow` | string | no | `mir-target` (default), `mir-lncrna`, or `mir-network` |
| `genome` | string | no | Species code (default `hg38`). Requesting `Targetscan` for a genome it has no dataset for returns HTTP 400 — see [Supported Genomes](#supported-genomes--species) |
| `cores` | int | no | CPU cores (default `ISOTAR_DEFAULT_CORES`, capped at `ISOTAR_MAX_CORES_PER_JOB`) |

**`mir-target` / `mir-lncrna`** — one miRNA per job:

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `mirna_id` | string | yes | miRNA identifier (e.g. `hsa-miR-21-5p`), or a free-form label when `mirna_seq` is given |
| `mirna_seq` | string | no | Raw custom sequence (17–30 nt, `ACGTU`). Bypasses the miRBase lookup; `modifications` / `shift` do not apply, and `mirna_id` is then restricted to `[A-Za-z0-9._-]` because it flows into output filenames |
| `modifications` | string[] | no | Nucleotide substitutions, `"<pos>:<from>\|<to>"` — e.g. `"8:A\|G"` (1-based) |
| `shift` | string | no | Re-cut from the precursor, `"left\|right"` — e.g. `"-2\|-3"` |
| `pre_id` | string | no | Precursor miRNA ID, for miRNAs with several precursors |
| `targets` | string[] | no | Restrict prediction to these targets. Gene symbols/accessions resolve to RefSeq; lncRNA IDs are normalised. Any unresolved target rejects the job with HTTP 400 |

**`mir-network`** — a list of miRNAs across both pools:

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `mirna_ids` | string[] | yes | Up to `ISOTAR_MAX_NETWORK_MIRNAS` miRNA IDs (deduped, order preserved) |
| `pairs` | object[] | no | ceRNA pairs `{gene, lncrna, corr}`. `corr` is optional and holds `mirna_gene`, `mirna_lncrna`, `gene_lncrna` — each optional, each a magnitude `0 <= x <= 1`. Omitted correlations default to `1.0` at scoring time and are listed in `corr_defaulted` on the response. Gene symbols are resolved to RefSeq up front |
| `restrict_to_pairs` | bool | no | Scan only the pairs' own targets instead of the whole reference (default `ISOTAR_NETWORK_RESTRICT_TO_PAIRS`). Requires `pairs`. Orders-of-magnitude faster, but the job no longer leaves behind a genome-wide table |
| `pre_ids` | object | no | `{"<mirna_id>": "<pre_id>"}` — only needed for multi-precursor miRNAs |
| `shifts` | object | no | `{"<mirna_id>": "left\|right"}` |
| `modifications` | object | no | `{"<mirna_id>": ["8:A\|U", ...]}` |
| `variants` | object | no | `{"<mirna_id>": [{"shift": "-7\|1"}, {"modifications": ["8:A\|U"]}, ...]}` — superset of the two above; each spec becomes its own graph node |

Keys in `pre_ids` / `shifts` / `modifications` / `variants` must name a miRNA in
`mirna_ids`, so a typo returns 400 instead of being silently dropped.

Response `202`:
```json
{ "job_id": "<uuid>", "task_id": "<celery-task-id>" }
```

---

### Get job status

```
GET /api/v1/jobs/<job_id>
```

Response includes `status` (`queued` | `running` | `succeeded` | `failed` |
`killed`) and a `progress` object showing per-tool status when the prediction
step is active. For `mir-network` jobs the two pools' progress files are merged
into one view. A job that succeeded with some tools failing carries
`tool_errors`.

---

### Query results

```
GET /api/v1/jobs/<job_id>/result
```

Returns 409 unless the job succeeded. The result database is built lazily on the
first call and cached as `result.db` in the output directory.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sortBy` | `tool_count` | Sort by `tool_count` or `gene_label` |
| `ascendOrDescend` | `desc` | `asc` or `desc` |
| `offset` | `0` | Pagination offset |
| `number` | `20` | Page size (max 1000) |
| `keyword` | — | Case-insensitive substring matched against gene label, gene name, **or** tool name |
| `pool` | — | `gene` or `lncrna`. **Required** for `mir-network` jobs, rejected for any other workflow |

Response:
```json
{
  "total_genes": 1500,
  "total": 42,
  "genes": [
    {
      "gene_id": "NM_001234",
      "gene_label": "TP53",
      "gene_name": "tumor protein p53",
      "tool_count": 4,
      "tools": ["DMISO", "PITA", "RNAhybrid", "Targetscan"]
    }
  ],
  "venn": {
    "sets": { "PITA": 310, "RNAhybrid": 220, "Targetscan": 180 },
    "intersections": { "PITA&RNAhybrid": 85, "PITA&Targetscan": 60 },
    "combinations": [{ "tools": ["PITA", "RNAhybrid"], "size": 85 }],
    "degrees": { "1": 420, "2": 145 }
  }
}
```

For a `mir-lncrna` job (or the `lncrna` pool of a network job) rows are keyed by
transcript ID with no gene-symbol mapping.

---

### Get the network graph

```
GET /api/v1/jobs/<job_id>/network
```

`mir-network` jobs only; 400 otherwise. Builds and caches `network.db` from both
pools, then shapes a tripartite graph.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `topGenes` | `50` | Discovery mode only — cap on top-degree genes |
| `topLncrna` | `50` | Discovery mode only — cap on top-degree lncRNAs |

With user pairs the response is in `pairs` mode and contains only the ceRNA
co-target bridges; without pairs it is a bounded `discovery` graph over the
top-degree nodes.

```json
{
  "mode": "pairs",
  "nodes": [{ "id": "NM_001234", "type": "gene", "label": "TP53", "name": "..." }],
  "edges": [{ "source": "NM_001234", "target": "hsa-miR-21-5p",
              "side": "gene", "tools": ["PITA"], "tool_count": 1 }],
  "pairs": [{ "gene": "NM_001234", "gene_input": "TP53", "gene_label": "TP53",
              "lncrna": "ENST00000...",
              "corr": { "mirna_gene": 0.9, "mirna_lncrna": 0.8, "gene_lncrna": 0.7 },
              "corr_defaulted": [], "corr_product": 0.504,
              "score": 0.130132, "score_version": 1,
              "best_bridge": { "mirna": "hsa-miR-21-5p,-7|1,shifted",
                               "score": 0.504, "is_variant": true },
              "bridges": [{ "mirna": "hsa-miR-21-5p", "is_variant": false,
                            "gene_tools": 2, "lncrna_tools": 1,
                            "consensus": 0.258199, "score": 0.130132 }],
              "bridge_mirnas": ["hsa-miR-21-5p"] }],
  "summary": { "mode": "pairs", "gene_count": 1, "mirna_count": 1,
               "lncrna_count": 1, "edge_count": 2, "pair_count": 1,
               "truncated": false, "score_version": 1,
               "tool_totals": { "gene": 6, "lncrna": 5 } }
}
```

#### ceRNA pair score (`score_version: 1`)

```
score = corr(miRNA,gene) x corr(miRNA,lncRNA) x corr(gene,lncRNA)
        x sqrt( (gene_tools / tool_totals.gene) x (lncrna_tools / tool_totals.lncrna) )
```

* **Correlations** are the caller's, as magnitudes in `0..1`. An omitted one
  contributes `1.0` and is named in `corr_defaulted`, so a pair with no
  correlation data scores on prediction consensus alone rather than scoring
  zero. Their product is reported separately as `corr_product`.
* **`tool_totals`** are the tools the job actually ran per pool — every selected
  tool on the gene side, minus TargetScan on the lncRNA side. Taken from the
  job's manifest, not hardcoded, so a species with no TargetScan release is
  scored out of 5 on the gene side instead of being capped below 1.
* **Scores are per bridging miRNA node**, in `bridges[]`, because tool counts
  live on edges and a miRNA's variants are separate nodes with their own
  support. Correlations are measured on the base miRNA and inherited unchanged
  by its variants, so within one pair `corr_product` is constant and WT-vs-variant
  differences come entirely from `consensus`.
* **The pair's `score` is the WT node's**, because a pair-level score is a claim
  about the real miRNA — an in-silico variant must not set it. It is `null` when
  WT does not bridge the pair at all; the pair is still returned, since "the
  variant bridges what WT cannot" is a result. `best_bridge` names the
  highest-scoring node either way.
* **`pairs` is sorted by `score`**, descending, with `null` scores last.

In `discovery` mode there are no user correlations, so all three default to
`1.0` and the score reduces to the consensus term. Discovery scores every
candidate before truncating to `MAX_DISCOVERY_PAIRS`, so the returned pairs are
the top-scoring ones rather than whichever the gene-major cartesian product
emitted first. Note that discovery ranks candidate targets by how many distinct
miRNA nodes hit them, which carries little signal for a single-miRNA job — the
isotar frontend therefore requires pairs and never reaches this mode.

---

### Enrichment analysis

```
POST /api/v1/jobs/<job_id>/enrichment
```

Runs `v2/enrichment_analysis.py` (gseapy/Enrichr) synchronously over a gene list
and writes its outputs into the job directory.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `genes` | string[] | yes | Gene symbols to test |
| `organism` | string | no | Default `Human` |
| `cutoff` | number | no | Adjusted p-value cutoff (default `0.05`) |

```
GET  /api/v1/jobs/<job_id>/enrichment           # parsed per-database reports
GET  /api/v1/jobs/<job_id>/enrichment/dotplot   # image/png
```

`GET .../enrichment` returns 404 until enrichment has been run, then
`{job_id, databases: {<gene_set>: [{term, overlap, p_value, adjusted_p_value,
odds_ratio, combined_score, genes}]}, has_dotplot}`.

---

### Literature verification (mir-network jobs)

```
GET  /api/v1/verification/config
POST /api/v1/jobs/<job_id>/verification
GET  /api/v1/jobs/<job_id>/verification
GET  /api/v1/jobs/<job_id>/verification/download   # text/tab-separated-values
```

Checks every predicted triplet of a succeeded `mir-network` job (one triplet per
ceRNA pair x bridge miRNA, ranked by pair then bridge score) against PubMed, and
optionally has a language model read the hits that could matter. Runs as a
Celery task on the dedicated `verification` queue (one process, so the whole
deployment stays under NCBI's 10 requests/second); progress and results live in
`<job>/verification/`. Human genomes only (`ISOTAR_VERIFICATION_GENOMES`).

`GET /verification/config` reports what the deployment can do:
`{enabled, llm_review_available, model, tissues, genomes, max_triplets, max_papers}`.
`enabled` is false until `NCBI_API_KEY` is set (POST then returns 503);
`llm_review_available` is false until an Anthropic-compatible gateway is
configured (POST with `llm_review: true` then returns 400).

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `tissue` | string | no | Tier 1 disease context: a preset (`lung` default, `breast`, `colorectal`, `gastric`, `liver`, `prostate`, `pancreatic`, `glioma`), `none` (tier 1 not run), or a raw PubMed boolean with `[tiab]`/`[mh]` fields |
| `llm_review` | bool | no | Also read the queued hits with the model (default `false`) |

POST returns 202 `{job_id, status: "queued", task_id, params}`; 409 while a
verification is already queued or running.

`GET .../verification` returns 404 until a verification has been started, then
the progress record `{job_id, status: queued|running|succeeded|failed, step,
arms_done, arms_total, triplets, queue_total, reviewed, error, params, ...}` with
`result` filled in once `status` is `succeeded`:

- `summary` -- `triplets`, `with_any_arm`, `with_two_arms`, `with_three_arms`,
  `same_paper_hits`, `unsearchable_lncrna`, `triplet_tier` counts, `truncated`.
- `rows` -- one per triplet, the published verification columns first:
  `id, miRNA, gene, lncRNA, triplet_tier, arms_supported, same_paper_hit,
  same_paper_url`, then `<arm>_pair, <arm>_tier, <arm>_url, <arm>_note` for
  `mir_gene`, `mir_lnc`, `gene_lnc`; then the site extras `gene_id, lncrna_id,
  mirna_node, pair_score, bridge_score, bridge_is_variant`.
- `queue` -- the hits handed to the model (same-paper hits, three-arm triplets,
  two-arm triplets with a tissue-tier arm or an experimental-looking title).
- `review` -- `{status: disabled|skipped|succeeded|failed, model, queue_size,
  reviewed, cached, errors, error, verdicts: [{id, scope, pmid, verdict,
  arm_match, tissue, evidence_type, sentences, note}]}`.

How to read a row: an arm tier is keyword co-occurrence of two symbols in one
title/abstract -- `1` with the tissue's cancer terms, `2` with any cancer term,
`3` symbols only, `none` no record, `unsearchable` when the lncRNA has no HGNC
symbol (its Ensembl id never appears in abstracts). A lower tier number records
disease context, not stronger evidence, and the three arms usually come from
three unrelated papers. `triplet_tier` is `tier{worst arm}` when all three arms
hit, `partial` with one or two, `none` otherwise. `same_paper_hit` is the only
field that points at a single study naming the whole axis, and only a verdict
of `reported_interaction` with its quoted sentence says the paper reports it.
The miRNA is searched at family level (arm suffix dropped), so an
`opposite_arm` verdict is common and worth showing.

lncRNA symbols come from `app_v1/lncrna_symbols.tsv`, built from the Ensembl
GRCh38 ncRNA dump by `reference_mapping/build_lncrna_symbols.py` (hg19 reuses
the symbols through the assembly-stable ENSG ids). The verdict step runs
`app_v1/verification/llm_review.py` under `python3.11` with the Anthropic SDK
pointed at the gateway in `ANTHROPIC_BASE_URL`; verdicts are cached by
(PMID, scope, symbols, model) in `ISOTAR_VERDICT_CACHE`.

---

### Download results

```
GET /api/v1/jobs/<job_id>/result/download
```

Returns a ZIP archive of all raw prediction output files, the run inputs
(`job.json`, `mirna.fa`, `targets.txt`), and — if enrichment analysis has been
run for the job — its outputs under `enrichment/` (`gene_list.csv`, one
`*.enrichr.reports.txt` per database, and `enrichment_dotplot.png`), and — if a
literature verification has completed — its tables under `verification/`
(`verification.tsv`, `verification.json`, `triplets.tsv`, `verdicts.tsv`).
Running or re-running either after a download invalidates the cached archive,
so the next download always reflects the latest results.

---

### Kill a job

```
POST /api/v1/jobs/<job_id>/kill
```

Terminates a queued or running job.

---

### Delete a job

```
DELETE /api/v1/jobs/<job_id>
```

Removes all files and metadata for the job.

---

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `ISOTAR_JOB_DIR` | `/opt/out/jobs` | Root directory for job storage |
| `ISOTAR_LOG_DIR` | `/opt/out` | Log file directory |
| `ISOTAR_LOG_LEVEL` | `INFO` | Logging level |
| `ISOTAR_VERSION` | (build arg) | Version reported by `/api/v1/version`; falls back to the bundled `VERSION` file |
| `ISOTAR_CORS_ORIGINS` | `*` | Comma-separated allowed origins for `/api/*` |
| `ISOTAR_RESOURCES_DIR` | `/opt/resources` | Per-species miRNA metadata JSON files |
| `ISOTAR_REFERENCE_MAPPING_DB` | `<app_v1>/reference_mapping.db` | Gene ID → symbol/name mapping database |
| `ISOTAR_LNCRNA_REFERENCE_DB` | `<app_v1>/lncrna_reference.db` | lncRNA transcript/gene ID reference database |
| `CELERY_BROKER_URL` | `amqp://` | RabbitMQ broker URL |
| `CELERY_RESULT_BACKEND` | `rpc://` | Celery result backend |
| `ISOTAR_MAX_CORES_PER_JOB` | `8` | Per-job core ceiling. The API rejects `cores` values above this with HTTP 400. |
| `ISOTAR_DEFAULT_CORES` | `8` | Default `cores` value when a job submission omits it. |
| `ISOTAR_MAX_CONCURRENT_JOBS` | `4` | Celery `--concurrency` for the `app_v1_worker` — caps the number of jobs running in parallel. |
| `ISOTAR_MAX_NETWORK_MIRNAS` | `20` | Max miRNAs in a single `mir-network` job |
| `ISOTAR_MAX_VALIDATE_TARGETS` | `5000` | Max targets per `POST /targets/validate` call |
| `ISOTAR_NETWORK_RESTRICT_TO_PAIRS` | `false` | Deployment default for a network job's `restrict_to_pairs` |
| `NCBI_API_KEY` / `NCBI_EMAIL` | (unset) | E-utilities credentials for literature verification. Verification is disabled until the key is set. Never logged. |
| `ISOTAR_VERIFICATION_GENOMES` | `hg19,hg38` | Genomes a verification may run on (the search terms and lncRNA symbol table are human) |
| `ISOTAR_VERIFICATION_QUEUE` | `verification` | Celery queue of the single-process verification worker (must match `app_supervisord.conf`) |
| `ISOTAR_VERIFICATION_DIR` | `<parent of ISOTAR_JOB_DIR>/verification` | Shared caches: `pubmed_cache.db` (E-utilities responses by query) and `verdicts.db` |
| `ISOTAR_PUBMED_CACHE` / `ISOTAR_VERDICT_CACHE` | (see above) | Override either cache file's path |
| `ISOTAR_PUBMED_CACHE_DAYS` | `0` | Refetch cached PubMed answers older than this many days (`0` = keep forever) |
| `ISOTAR_PUBMED_RETMAX` | `20` | PMIDs kept per query |
| `ISOTAR_VERIFICATION_MAX_TRIPLETS` | `500` | Top-ranked triplets searched per job (`summary.truncated` says when the cap applied) |
| `ANTHROPIC_BASE_URL` / `ANTHROPIC_AUTH_TOKEN` | (unset) | Anthropic-compatible gateway (e.g. LiteLLM) and its key, read by the SDK in `llm_review.py`. Model review is offered only when the token is set. |
| `ISOTAR_VERIFICATION_MODEL` | `claude-opus-5-5` | Model name as the gateway exposes it |
| `ISOTAR_VERIFICATION_MAX_PAPERS` | `40` | Papers read by the model per job (cached verdicts are free and not counted) |
| `ISOTAR_VERIFICATION_TIMEOUT` | `3600` | Seconds allowed for the model review step |
| `ISOTAR_VERIFICATION_PYTHON` | `python3.11` | Interpreter for `llm_review.py` (needs the `anthropic` and `pydantic` packages) |
| `ISOTAR_LNCRNA_SYMBOLS_TSV` | `<app_v1>/lncrna_symbols.tsv` | lncRNA gene id -> HGNC symbol table |

## Project Structure

```
isoTar-v2/
├── app/                        # Legacy Flask application (Python 2.7)
├── app_v1/                     # REST API v1 (Python 3.6)
│   ├── app.py                  # Flask routes and the run_job Celery task
│   ├── celery_app.py           # Celery configuration
│   ├── limits.py               # Core / concurrency limits and validation
│   ├── version.py              # Version resolution for /api/v1/version
│   ├── result_db.py            # Gene result database (build + query + Venn)
│   ├── lncrna_results.py       # Same, for lncRNA-pool results
│   ├── network_results.py      # Tripartite network graph builder
│   ├── target_resolver.py      # Gene symbol / accession → RefSeq resolution
│   ├── lncrna_reference.py     # lncRNA ID validation and normalisation
│   ├── parse_result.py         # Tool output parsers
│   ├── logger.py               # Date-rotating file logger
│   ├── verification/           # Literature verification (PubMed tiers, roll-up, model review)
│   │   ├── pubmed.py           # Query terms, rate-limited client, SQLite cache
│   │   ├── rollup.py           # Per-arm tiers -> per-triplet rows
│   │   ├── review.py           # Which hits the model reads
│   │   ├── symbols.py          # lncRNA transcript/gene id -> HGNC symbol
│   │   ├── papers.py           # Abstract / full-text fetch, sentence extraction
│   │   ├── pipeline.py         # Celery-side driver, writes <job>/verification/
│   │   └── llm_review.py       # python3.11: verdict per paper via the Anthropic SDK
│   ├── lncrna_symbols.tsv      # lncRNA gene id -> symbol (reference_mapping/build_lncrna_symbols.py)
│   ├── reference_mapping.db    # Gene ID → label/name mapping
│   └── lncrna_reference.db     # lncRNA transcript/gene ID reference
├── v2/                         # Prediction engine scripts (run as subprocesses)
│   ├── mirna_processing.py     # miRNA FASTA generation (python2.7)
│   ├── mirna_predicting.py     # Tool orchestration and parallelism
│   ├── parse_result.py         # Result parsing library
│   ├── lncrna_ids.py           # lncRNA header/ID helpers
│   ├── targetscan_new.py       # Standalone TargetScan runner
│   ├── process_predictions.py  # Ad-hoc prediction post-processing
│   ├── enrichment_analysis.py  # Gene enrichment analysis (gseapy/Enrichr)
│   └── tests/                  # unittest suite for the v2 scripts
├── v2/opt/reference_files/     # 3' UTR FASTA files per species
├── v2/opt/resources/           # Per-species miRNA metadata JSON
├── reference_mapping/          # Sources and builders for the reference DBs
├── scripts/                    # release.sh, reference builders, helpers
├── Dockerfile                  # Production image (builds on base)
├── isotar-base.Dockerfile      # Base image (Ubuntu 16.04 + Python 3.6 / 3.11)
├── .env.template               # Every env var the server reads (copy to .env, pass with --env-file)
├── app_supervisord.conf        # Supervisor process config
├── nginx.conf / app.conf       # Nginx configuration
└── CLAUDE.md / AGENTS.md       # Developer notes (kept identical)
```

## Development

### Tests

The `v2/tests/` suite is plain `unittest`. There is no `__init__.py`, so it runs
as an implicit namespace package from the repo root — name the modules
explicitly rather than using `unittest discover`:

```bash
python3 -m unittest v2.tests.test_parse_result          # one module
python3 -m unittest $(ls v2/tests/test_*.py | sed 's|/|.|g; s|\.py$||')   # all
```

Modules that import `app_v1` need Flask and Celery installed; without them they
error out at import. Everything else runs against the standard library alone.

### Git pre-commit hook (one-time per clone)

This repo ships a pre-commit hook that blocks commits introducing non-ASCII
bytes into `v2/*.py` without a PEP 263 encoding declaration (Python 2.7
otherwise fails with `SyntaxError: Non-ASCII character ...`). Enable it once
after cloning:

```bash
git config core.hooksPath .githooks
```

The hook only runs when `v2/*.py` files are staged, so it's silent on
unrelated commits. Bypass in emergencies with `git commit --no-verify`.


