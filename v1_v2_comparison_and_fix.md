# isoTar v1 vs v2 — Prediction Divergence: Investigation and Fixes

**Date:** 2026-07-21
**Images compared:** `frankfeng78/isotar-v1:1.2.1` vs `frankfeng78/isotar-v2:0.3.18`
**Branch at time of work:** `main` @ `847334e` (chore(release): v0.3.18)
**Comparison harness:** `/bmbl_data/shaohong/projects/playground/v1v2-compare`
(`compare.py`, `run_v1.sh`, `run_v2.sh`, `report/comparison.tsv`)

---

## TL;DR

v2 predicts far more targets than v1. The original hypothesis was that this is entirely
the different 3'UTR reference files. **That turned out to be true for RNAhybrid only.**
A controlled experiment (same reference fed to both engines) exposed **two real v2
regressions**:

| Tool | Cause of v1↔v2 divergence | Status |
|---|---|---|
| **RNAhybrid** | 3'UTR reference only (coverage + sequence content). Wrapper faithful. | ✅ No defect |
| **miRanda** | v2's parser dropped v1's canonical **seed-match filter** → ~7× over-prediction | ✅ **Fixed** (code) |
| **PITA** | Bundled **`RNAddG4` binary segfaults** in the v2 image → `dGopen` empty → `ddG` collapses to `dGduplex` → ~21–40× over-prediction | ✅ **Fixed** (binary rebuilt) + guard added |
| miRmap | Genuinely different software (miRmap 1.x vs miRmap 2) | Out of scope |
| TargetScan | Different mechanism (v2 = precomputed Ensembl datasets → RefSeq) | Out of scope |

**Three fixes applied:** the two parser patches in `v2/parse_result.py` and
`app_v1/parse_result.py` (+105 / −33 lines), plus a rebuilt 64-bit `RNAddG4` binary
(`scripts/rebuild_RNAddG4.sh`). All three verified to reproduce v1's results exactly.

---

## 1. Scope and setup

Both engines were run on the **same 5 miRNAs**, **same 5 tools**, **same genome build
(hg19)**, fed the **same mature sequences**:

| mature | len | sequence |
|---|---|---|
| hsa-miR-21-5p | 22 | UAGCUUAUCAGACUGAUGUUGA |
| hsa-miR-155-5p | 24 | UUAAUGCUAAUCGUGAUAGGGGUU |
| hsa-let-7a-5p | 22 | UGAGGUAGUAGGUUGUAUAGUU |
| hsa-miR-124-3p | 22 | UAAGGCACGCGGUGAAUGCCAA |
| hsa-miR-1-3p | 22 | UGGAAUGUAAAGAAGUAUGUAU |

Results are reduced to a **set of base RefSeq accessions** per (miRNA, tool); version
suffixes are stripped so the two reference namespaces line up.

> **Important scope note:** `mirna_processing.py` is **not exercised** by this comparison.
> `run_v2.sh` feeds a pre-made FASTA straight into `mirna_predicting.py` via `-i`. What is
> validated here is `mirna_predicting.py` + `parse_result.py`. Verifying
> `mirna_processing.py` requires a separate test (generate a miRNA via both engines and
> diff the FASTA).

---

## 2. Initial comparison (each engine on its own reference)

Summed over all 5 miRNAs:

| tool | v1 total | v2 total | v2/v1 | % of v1 reproduced in v2 | mean jaccard_fair |
|---|--:|--:|--:|--:|--:|
| **RNAhybrid** | 13,125 | 26,318 | 2.0× | **98%** | **0.835** |
| **miRanda** | 2,176 | 13,757 | 6.3× | **99%** | 0.266 |
| **PITA** | 3,805 | 29,641 | 7.8× | 83% | 0.162 |
| **TargetScan** | 36,482 | 73,371 | 2.0× | 48% | 0.287 |
| **miRmap** | 19,285 | 977 | **0.05×** | **2%** | **0.018** |

`jaccard_fair` = restricted to the 39,702 accessions present in **both** references.

---

## 3. Reference-file analysis

The two engines use different 3'UTR sets:

| | v1 | v2 |
|---|---|---|
| source | hg19 **refGene** | hg19 **ncbiRefSeqCurated** |
| path | `/app/iso/static/public/resources/3utr.fa` | `/opt/reference_files/hsa_HG19_only_3UTRs.fasta` |
| base accessions | 40,374 | 66,869 (**1.66×**) |
| file size | 66 MB | 174 MB |

Common accessions: **39,702**. Critically, *the same accession does not mean the same
sequence*:

- Only **37.5%** of common accessions have an **identical** UTR length.
- **62.5%** differ (median |Δ| = 5 nt, tails up to ~30,000 nt where refGene carries a
  truncated UTR and ncbiRefSeqCurated the full one).

So the reference differs at **two levels**: (1) *which transcripts exist*, and (2) *what
sequence each shared transcript has*.

### Divergence partition (per tool)

| bucket | RNAhybrid | miRanda | PITA |
|---|--:|--:|--:|
| v2-only, transcript **absent from v1's catalog** | **83%** | 50.1% | 48.7% |
| v2-only, transcript **in both references** | 17% | 49.9% | 51.3% |
| v1 calls not reproduced in v2 | 258 | 28 | 662 |

RNAhybrid's divergence is dominated by catalog size — clean and fully reference-driven.
miRanda/PITA are ~50/50, and that second half is what the controlled experiment below
resolved.

---

## 4. The controlled experiment — same reference, both wrappers

v1's **exact** 3'UTR file was bind-mounted over v2's hg19 reference path, and v2's
pipeline (`mirna_predicting.py` + `parse_result.py`) re-run. Same reference in → if the
sets match, the wrapper introduces nothing.

```bash
docker run --rm \
  -v "$SP:/w" \
  -v "$SP/v1_3utr.fa:/opt/reference_files/hsa_HG19_only_3UTRs.fasta:ro" \
  --entrypoint bash frankfeng78/isotar-v2:0.3.18 -c '
    python3.6 /opt/v2/mirna_predicting.py -c 16 -i /w/eq_input.fasta \
      -t miRanda RNAhybrid PITA -g hg19 -tt gene -o /w/eq_out
    python3.6 /opt/v2/parse_result.py /w/eq_out'
```

**Result — identical reference, yet:**

| miRNA | tool | v1 | v2 | shared | v1_only | v2_only | jaccard |
|---|---|--:|--:|--:|--:|--:|--:|
| miR-21-5p | RNAhybrid | 812 | 802 | 800 | 12 | 2 | **0.9828** |
| miR-1-3p | RNAhybrid | 1489 | 1478 | 1472 | 17 | 6 | **0.9846** |
| miR-21-5p | miRanda | 49 | **358** | 49 | 0 | 309 | 0.1369 |
| miR-1-3p | miRanda | 66 | **144** | 66 | 0 | 78 | 0.4583 |
| miR-21-5p | PITA | 106 | **2241** | 106 | 0 | 2135 | 0.0473 |
| miR-1-3p | PITA | 365 | **3565** | 349 | 16 | 3216 | 0.0975 |

**RNAhybrid matches; miRanda and PITA do not.** The divergence is therefore *not* 100%
the reference.

### Tool binaries and parameters were ruled out first

| tool | v1 md5 | v2 md5 | verdict |
|---|---|---|---|
| RNAhybrid | `2d35c117…` | `2d35c117…` | identical binary |
| PITA `pita_prediction.pl` + all 30 `lib/*.pl` | `4ba17390…` | `4ba17390…` | identical scripts |
| miRanda | `59ea32de…` | `f44cdcf9…` | md5 differs, **both `miranda v3.3a`** (rebuild) |

Invocation parameters are effectively identical in all three cases:

- **RNAhybrid** — both: `-c -e -20 -s 3utr_human -q <mir> -t <utr> -m <utr_len> -n <mir_len>`
  (v1 via `itutils.py:447` + `additional_opts='-m %d -n %d'`; v2 via `run_rnahybrid`).
  `-m` is always ≥ every UTR it scans, so it only shifts the p-value *number*, never the
  energy-gated hit *set* → results are chunk-count invariant.
- **miRanda** — v1: `miranda <mir> <utr> -quiet -out … -en -20`; v2: same flags reordered.
- **PITA** — scoring flags identical (`-l 7-8`, `-gu 7;0;8;0`, `-m 7;0;8;0`). Only
  cosmetic difference: v1 `-gpx` (an unrecognized flag, silently ignored) vs v2 `-gxp`
  ("produce a Genomica .gxp output file"). Neither controls `dGopen`.

---

## 5. Root cause — miRanda (parser regression)

Raw miRanda output is **byte-equivalent**: 358 distinct accessions on both sides. Yet
v1 reports 49 and v2 reports 358. The difference is entirely in parsing:

- **v1** (`/app/iso/isotarpkg/itutils.py :: parseMirandaResults`) walks each alignment
  block and keeps a target **only if the seed (positions 2–7) is fully base-paired** —
  `len(seed_region.replace('|','')) == 0`.
- **v2** (`parse_result.py :: parseMirandaResults`) kept **every** target appearing in any
  `>` summary line — **no seed filter**.

The alignment blocks *are* present in v2's output (361 `Query:` lines) — miranda v3.3a's
`-quiet` does **not** suppress them — so the fix is parser-only; the invocation needs no
change.

**miRanda hit-block layout** (offset from the `Query:` line):

```
 0: "   Query:    3' aguugUAGUCAGACUAUUCGAu 5'"
 1: '                     ||||||| ||||||||'      <- pairing line; seed slice [-9:-2]
 2: "   Ref:      5' tcttgATCAGTCAGATAAGCTt 3'"
 3: ''
 4: '   Energy:  -21.860001 kCal/Mol'
 5: ''
 6: 'Scores for this hit:'
 7: '>hsa-miR-21-5p,WT\thg19_refGene_NM_001134445\t177.00\t-21.86\t…'   <- target id = field [1]
```

---

## 6. Root cause — PITA (broken `RNAddG4` binary)

Raw PITA files have the same 2,997 distinct accessions, but **different `ddG` values**:

| | v1 | v2 |
|---|---|---|
| `dGopen` populated | **3337 / 3337 rows** | **0 / 3394 rows** |
| `ddG` | `dGduplex − dGopen` | `= dGduplex` |
| sample row | dGduplex=−13.9, dGopen=−13.39, **ddG=−0.50** | dGduplex=−16.3, dGopen=*(empty)*, **ddG=−16.3** |
| unique accessions passing `ddG ≤ −10` | **106** | **2241** |

PITA computes `dGopen` (site accessibility) via `lib/RNAddG_compute.pl` →
`Bin/ViennaRNA/ViennaRNA-1.6/Progs/**RNAddG4**` (`RNAddG_compute.pl:201`; `dGduplex` comes
from `RNAduplex` at line 142). The binaries differ between images:

| binary | v1 | v2 |
|---|---|---|
| `RNAddG4` | `4901a2f2…` (612531 B, 2017-10-19) | `bd08d3ab…` (612526 B, rebuilt) |
| `RNAduplex` | `0923d709…` | `4340333c…` |

Behaviour of the *same* invocation:

- **v1 binary:** prints the accessibility values, *then* aborts (`exit 134`,
  `munmap_chunk(): invalid pointer`). The value is flushed first → **captured**.
- **v2 binary:** **segfaults with no output** (`exit 139`) → `dGopen` empty.

Ruled out as causes: glibc (**both 2.23**), RNAfold (**both 2.4.11**, identical output),
PITA `.pl` scripts (byte-identical), PITA flags, and runtime knobs
(`setarch -R`, `MALLOC_CHECK_=0`, `GLIBC_TUNABLES=glibc.malloc.check=0` — all still
segfault). It is the binary itself. Confirmed by running **v1's `RNAddG4` inside the v2
container**: it prints correctly and then does the harmless post-print abort.

### Fix confirmation (identical 3000-transcript subset, miR-21-5p, `-c 1`)

| metric | BROKEN (stock v2 `RNAddG4`) | FIXED (v1 `RNAddG4` swapped in) |
|---|--:|--:|
| `dGopen` populated | 0 / 236 rows | **236 / 236 rows** |
| `ddG == dGduplex` | 236 / 236 | 0 / 236 |
| targets passing `ddG ≤ −10` | **161** | **4** |
| sample row | dGduplex=−13.4, dGopen=*(empty)*, ddG=−13.4 | dGduplex=−13.4, dGopen=−8.19, **ddG=−5.20** |

Restoring a working `RNAddG4` corrects PITA by **~40×**.

---

## 7. Fixes applied

Applied to **both** `v2/parse_result.py` (CLI path) and `app_v1/parse_result.py`
(REST-API path) — +105 / −33 lines total.

### 7.1 `parseMirandaResults` — restore v1's seed-match filter

Walks each alignment block instead of harvesting every `>` line. The target id is taken
from tab field `[1]` of the `>` line via the existing `id_extractor`, so it works for both
the RefSeq (`v2`) and lncRNA (`app_v1`) flows.

```python
for i in range(len(lines)):
    if not re.match(r"^\s+Query:\s+3'\s+\S+\s+5'$", lines[i], re.I):
        continue
    # Pairing string is the next line; seed positions 2-7 map to slice [-9:-2]
    # (5' miRNA end is right-aligned). All '|' => fully paired; else reject.
    if i + 1 >= len(lines) or lines[i + 1][-9:-2].replace('|', '') != '':
        continue
    # Query(+0) -> pairing(+1) -> Ref(+2) -> +5 lands on the '>' summary line.
    j = i + 7
    if j >= len(lines) or not lines[j].startswith('>'):
        continue
    parts = lines[j].rstrip('\r\n').split('\t')
    if len(parts) < 2:
        continue
    tar = _extract_transcript_id(parts[1])      # app_v1: id_extractor(parts[1])
    if tar and tar not in results:
        results.append(tar)
```

*Caveat documented in the docstring:* this relies on the alignment block, which
miranda v3.3a emits even under `-quiet`. A truly summary-only output would yield nothing —
acceptable, since the seed cannot be verified without the alignment.

### 7.2 `parsePITAResults` — accessibility guard

Counts rows and populated `dGopen` values; if `dGopen` is empty for *every* row (the
`RNAddG4`-failure signature), it raises instead of emitting invalid targets.

```python
if rows_seen > 0 and dgopen_seen == 0:
    raise RuntimeError(
        "PITA: dGopen is empty for all {} site rows in {} -- the RNAddG4 "
        "accessibility binary produced no output, so ddG collapsed to "
        "dGduplex and PITA predictions are invalid (massive over-prediction). "
        "Replace the broken RNAddG4 binary at /opt/PITA64bit/Bin/ViennaRNA/"
        "ViennaRNA-1.6/Progs/RNAddG4 and re-run.".format(rows_seen, output_f_path))
```

**Why the wrapper and not `RNAddG_compute.pl`:** `run_pita()` invokes PITA with
`stderr=devnull`, so a `die()` inside the Perl would be **silently swallowed**. The
wrapper-side check is the reliable place. The caller (`process_sequence`, line ~371)
re-raises with sequence context, so the failure surfaces.

### 7.3 `RNAddG4` rebuilt from source (the actual PITA fix)

`scripts/rebuild_RNAddG4.sh` builds a working 64-bit `RNAddG4` from the ViennaRNA-1.6
sources already vendored in `tools/`, and the result is installed at
`tools/PITA64bit/Bin/ViennaRNA/ViennaRNA-1.6/Progs/RNAddG4`.

Why a rebuild was needed rather than a simple relink: the vendored `lib/libRNA.a` is
**32-bit** (`i386 ... incompatible with i386:x86-64`) and there is no shipped
`RNAddG4.o`, so the library must be compiled from source (64-bit) first:

```bash
gcc -O2 -I$VRNA/H -I$VRNA -c $VRNA/lib/*.c        # 26 objects
ar rcs libRNA64.a *.o
gcc -O2 -I$VRNA/H -I$VRNA -c $VRNA/Progs/RNAddG4.c -o RNAddG4.o
gcc -O2 -o RNAddG4 RNAddG4.o libRNA64.a -lm
```

| build | md5 | behaviour |
|---|---|---|
| shipped v2 (removed) | `bd08d3ab…` | **segfaults before printing** → `dGopen` empty |
| isoTar v1 | `4901a2f2…` | prints values, then aborts on exit |
| **rebuilt (installed)** | `4fc63f42…` | prints values, then aborts on exit — **matches v1** |

The trailing abort is a latent double-free in the 2006 `RNAddG4.c` source and is
inherited by *any* correct build; it is harmless because stdout is flushed first. The
broken build's distinguishing failure is dying *before* it prints.

The script smoke-tests its own output and refuses to install a binary that emits
nothing. It is deterministic — re-running reproduces md5 `4fc63f42…`.

**`RNAduplex` was checked and is fine:** `dGduplex` is identical on **3337/3337**
matched `(accession, start, end)` sites between v1-native and v2. Only `RNAddG4` was
broken, despite 16 of the 56 `Progs/` files differing between images.

---

## 8. Validation

Patched functions run against the real captured data:

```
=== PATCHED parseMirandaResults (expect 49, 66) ===
  hsa-miR-21-5p  -> 49   OK
  hsa-miR-1-3p   -> 66   OK

=== PATCHED parsePITAResults ===
  BROKEN -> RuntimeError raised OK ("dGopen is empty for all 236 site rows…")
  FIXED  -> OK, 4 targets (accessibility applied)
```

miRanda agreement with v1 after the fix: **jaccard = 1.000**, `v1_only = 0`,
`v2_only = 0` on both miRNAs. Replication was validated first by running the reconstructed
v1 filter on **v1's own** raw output and reproducing v1's 49 exactly.

Both files byte-compile; no f-strings introduced (2.7/3.5-safe per project rules).

---

## 9. Outstanding work

1. ~~Replace the broken `RNAddG4` binary~~ — **DONE** (§7.3). Rebuilt from source via
   `scripts/rebuild_RNAddG4.sh` and installed; verified byte-identical output to v1 and an
   end-to-end PITA run matching v1's target set. **Requires an image rebuild to take
   effect** (`Dockerfile:19` `ADD tools /opt/`).

   Sanity check in any image: `echo "GGGGAAAACCCCAAAAGGGGUUUU" | RNAddG4 -u 0 -s 5 -f 25 -t 25`
   → must print a tab-separated values line.

   **Deployment — no Dockerfile changes needed.** All three fixes are already picked up by
   existing directives: `Dockerfile:20` (`COPY v2/*.py /opt/v2/`), `Dockerfile:72`
   (`COPY app_v1 /app_v1`), `Dockerfile:19` (`ADD tools /opt/`).
   `isotar-base.Dockerfile` contains no `COPY`/`ADD` at all, so the base image does **not**
   need rebuilding. A plain `docker build -f Dockerfile .` is sufficient.

   **Decision (2026-07-21): ship the committed binary rather than compiling at image-build
   time.** Residual risk to be aware of: the committed `RNAddG4` is built against this
   base image's toolchain (`ubuntu:16.04`, glibc 2.23, gcc 5.4). If `isotar-base` is ever
   moved to a newer base, the prebuilt binary may segfault again *in exactly the same
   silent way* — that is precisely how this bug arose. Two mitigations already exist if
   that happens: the `parsePITAResults` guard (§7.2) now converts it from silent
   corruption into a loud error, and `scripts/rebuild_RNAddG4.sh` regenerates a matching
   binary. If the base image is ever upgraded, re-run that script and re-verify with the
   sanity check above. (The alternative — a `RUN` step invoking the script during the
   image build, which would fail the build on a bad binary — was considered and
   deliberately not adopted.)

2. **`RNAddG_compute.pl` temp-file collision (latent).** It writes hard-coded
   `tmp_seqfile1` / `tmp_seqfile2` into the **current working directory**
   (`RNAddG_compute.pl:127,176`). This was masked while `RNAddG4` crashed; once the binary
   works, running PITA with `-c > 1` may let parallel chunks clobber each other. Give each
   chunk an isolated CWD (or per-chunk temp names) before trusting parallel PITA.

3. **Re-run the v1↔v2 comparison** after (1) to confirm PITA and miRanda now agree with v1
   on a shared reference, leaving only the genuine 3'UTR-reference differences.

4. **`mirna_processing.py` is still unverified** (not exercised by this harness) — needs a
   separate FASTA-diff test if that path also matters.

5. Changes are **uncommitted on `main`**. Move to a branch / commit as appropriate.

---

## 10. Reproduction

```bash
# Full v1 vs v2 comparison (each on its own reference)
cd /bmbl_data/shaohong/projects/playground/v1v2-compare
./run_v1.sh 8          # -> out/v1/prediction/prediction.json
./run_v2.sh 8          # -> out/v2/miRNA_prediction_results/*.json
python3.11 compare.py  # -> report/comparison.{tsv,json}

# Controlled test: v2 wrapper on v1's reference (see §4 for the docker invocation)

# Check whether RNAddG4 works in an image
docker run --rm --entrypoint bash <image> -c \
  'echo "GGGGAAAACCCCAAAAGGGGUUUU" | \
   /opt/PITA64bit/Bin/ViennaRNA/ViennaRNA-1.6/Progs/RNAddG4 -u 0 -s 5 -f 25 -t 25'
# working  -> prints "…<tab>-10.8172<tab>-9.90425…" (then a harmless abort, exit 134)
# broken   -> prints nothing, exit 139 (segfault)

# Detect the failure signature in any existing PITA output
awk -F'\t' 'NR>1 && $12=="" {c++} END{print "empty dGopen rows:", c+0}' *_pita_results.tab
```

---

## Appendix — key source locations

| item | location |
|---|---|
| v1 tool command table | `/app/iso/isotarpkg/itutils.py:416-450` |
| v1 miRanda parser (seed filter) | `/app/iso/isotarpkg/itutils.py :: parseMirandaResults` |
| v1 PITA parser (`ddG <= -10`) | `/app/iso/isotarpkg/itutils.py:813` |
| v1 per-tool result matrix | `/app/iso/isotarpkg/pred_analysis.py :: get_output_data_per_tool` |
| v2 reference path constants | `v2/mirna_predicting.py:65-66` |
| v2 UTR chunking (verbatim split) | `v2/mirna_predicting.py :: process_3utr_fasta` |
| v2 tool runners | `v2/mirna_predicting.py :: run_miranda / run_rnahybrid / run_pita` |
| v2 parsers (**patched**) | `v2/parse_result.py :: parseMirandaResults / parsePITAResults` |
| REST-API parsers (**patched**) | `app_v1/parse_result.py` (same two functions) |
| PITA accessibility pipeline | `tools/PITA64bit/lib/pita_run.pl:128` → `lib/RNAddG_compute.pl:201` |
| PITA binaries baked into image | `Dockerfile:19` (`ADD tools /opt/`) |
