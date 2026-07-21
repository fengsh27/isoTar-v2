# isoTar v1 vs v2

## Why the two versions disagree — and what we found

`frankfeng78/isotar-v1:1.2.1`  vs  `frankfeng78/isotar-v2:0.3.18`

2026-07-21

---

## TL;DR

We assumed v1↔v2 differences came from the **different 3'UTR reference files**.

**That is true for only one of the five tools.**

| Tool | v2/v1 | Real cause | Verdict |
|---|--:|---|---|
| **RNAhybrid** | 2.0× | 3'UTR reference | ✅ No defect |
| **miRanda** | 6.3× | v2 parser dropped the seed filter | ❌ **Fixed** |
| **PITA** | 7.8× | Broken bundled binary | ❌ **Fixed** |
| **TargetScan** | 2.0× | v2 parser maps ENST→RefSeq by gene symbol | ❌ **Open** |
| **miRmap** | **0.05×** | v2 parser applies a ΔG gate v1 never had | ❌ **Open** |

**4 of 5 tools carry real v2 defects.** Three over-predict; miRmap under-predicts ~99%.

---

## The pattern

In **every** case, the prediction tool itself ran correctly.

- miRanda — raw output **identical** between versions
- TargetScan — raw output **byte-identical** (same perl, same dataset)
- miRmap — at engine level v2 loses **nothing** (`v1_only = 0`)
- PITA — the `.pl` scripts are byte-identical; a vendored helper **binary** was broken

> **Every defect is in isoTar's own integration layer** — `parse_result.py`, the wrapper,
> or the packaged binaries — **not in the science.**

Three of the four are post-processing bugs. None produced an error message.

---

## What we compared

- **Same 5 miRNAs** — miR-21-5p, miR-155-5p, let-7a-5p, miR-124-3p, miR-1-3p
- **Same 5 tools**, **same genome build (hg19)**, **same mature sequences**

Results reduced to a **set of RefSeq accessions** per (miRNA, tool), versions stripped.

Each suspicious tool then got a **controlled experiment**: hold the reference fixed,
run both engines, and see whether the divergence survives.

> `mirna_processing.py` is **not** exercised — the harness feeds a pre-made FASTA
> directly to `mirna_predicting.py`. This validates prediction + parsing only.

---

## The raw observation

Summed over all 5 miRNAs:

| Tool | v1 | v2 | v2/v1 | % of v1 reproduced |
|---|--:|--:|--:|--:|
| RNAhybrid | 13,125 | 26,318 | 2.0× | 98% |
| miRanda | 2,176 | 13,757 | **6.3×** | 99% |
| PITA | 3,805 | 29,641 | **7.8×** | 83% |
| TargetScan | 36,482 | 73,371 | 2.0× | **48%** |
| miRmap | 19,285 | **977** | **0.05×** | **2%** |

Four tools disagree wildly, in **both** directions. That alone rules out a single
shared cause.

---

## The references really are different

| | v1 | v2 |
|---|---|---|
| source | hg19 **refGene** | hg19 **ncbiRefSeqCurated** |
| transcripts | 40,374 | 66,869 (**1.66×**) |

And for the **39,702 accessions present in both**:

- Only **37.5%** have an identical UTR sequence length
- **62.5%** differ — same accession, *different* annotated 3'UTR

A real effect, at two levels: *which transcripts exist*, and *what sequence each has*.

**But it explains only RNAhybrid.**

---

## The key experiment

**Hold the reference fixed. Run both engines. Does the divergence survive?**

```bash
docker run --rm \
  -v "$SP/v1_3utr.fa:/opt/reference_files/hsa_HG19_only_3UTRs.fasta:ro" \
  frankfeng78/isotar-v2:0.3.18 \
  python3.6 /opt/v2/mirna_predicting.py -t miRanda RNAhybrid PITA -g hg19 ...
```

This single move separated "reference difference" from "software defect" —
and it is what turned a vague 2–8× discrepancy into four specific bugs.

---

## Result: the tools split

**Identical reference, both engines:**

| miRNA | Tool | v1 | v2 | Jaccard |
|---|---|--:|--:|--:|
| miR-21-5p | **RNAhybrid** | 812 | 802 | **0.98** ✅ |
| miR-1-3p | **RNAhybrid** | 1,489 | 1,478 | **0.98** ✅ |
| miR-21-5p | miRanda | 49 | **358** | 0.14 ❌ |
| miR-1-3p | miRanda | 66 | **144** | 0.46 ❌ |
| miR-21-5p | PITA | 106 | **2,241** | 0.05 ❌ |
| miR-1-3p | PITA | 365 | **3,565** | 0.10 ❌ |

TargetScan and miRmap were confirmed the same way in their own studies.

---

## Cause 1 — RNAhybrid: genuinely the reference ✅

No defect. Binary **md5-identical** between images; parameters identical
(`-c -e -20 -s 3utr_human -m … -n …`); 98% of v1's calls reproduced.

Divergence decomposes as:

- **83%** — transcripts absent from v1's smaller catalog
- **17%** — shared transcripts whose UTR *sequence* differs

**Action: none.** Expected behaviour from a newer reference.

*This is our control: it proves the harness and the reference hypothesis both work —
which is exactly why the other four results are damning.*

---

## Cause 2 — miRanda: v2 dropped the seed filter ❌

Raw miRanda output is **identical** (358 accessions both sides). The difference is parsing.

| | Behaviour | Result |
|---|---|--:|
| **v1** | keeps a target only if the miRNA **seed is fully Watson-Crick paired** | 49 |
| **v2 (shipped)** | keeps **every** `>` summary line — no seed check | **358** |

v1 inspects the alignment block; v2 harvested summary lines only. The seed-match test
was not carried over when `parse_result.py` was rewritten.

**Impact: ~7× over-prediction**, independent of the reference.
Includes sites with mismatched or G:U-wobbled seeds.

---

## Cause 3 — PITA: a broken bundled binary ❌

PITA scores **ΔΔG = ΔGduplex − ΔGopen**. ΔGopen (site accessibility) *is* PITA's premise.

1. v2 ships a **different build** of `RNAddG4` (computes ΔGopen) that **segfaults before printing**
2. `RNAddG_compute.pl` **never checks** exit status or output
3. It computes `ddG = dGduplex + (dG1 - dG0)`; undefined → Perl treats as `0`
4. → **`ddG = dGduplex`** — accessibility silently dropped
5. → rows clearing `ddG ≤ −10` jump from **3.2% → ~75%**

**PITA stopped being PITA** — it degenerated into "anything with good duplex energy."

---

## Cause 3 — PITA: the numbers

hsa-miR-21-5p:

| Stage | Targets | Cause |
|---|--:|---|
| v1 baseline | 106 | — |
| v2, **same** reference | **2,241** | broken binary (**21×**) |
| v2, own reference | 4,357 | + reference (~2×) |

**Ruled out** — all verified, not assumed:

- PITA `.pl` scripts — byte-identical (all 31 files)
- PITA parameters and the `ddG ≤ −10` threshold — identical
- `RNAduplex` — fine (dGduplex identical on **3337/3337** sites)
- glibc (2.23), RNAfold (2.4.11) — identical

---

## Cause 4 — TargetScan: a gene-symbol join ❌

**TargetScan never reads either 3'UTR reference.** Both engines run the same perl on the
same 6.8 GB pre-aligned dataset with the same miRNA input → **byte-identical raw output.**

**100% of divergence is in `parse_result.py`.** Two differences, opposite directions:

| | Difference | Direction | Magnitude |
|---|---|---|---|
| **A** | v2 filters rows to `species_ID == 9606`; v1 does not | v2 **smaller** | 14,124 → 8,426 ENST (−40%) |
| **B** | v1 maps ENST→RefSeq **1:1** (curated table); v2 joins **on gene symbol**, 1:many | v2 **larger** | mean fanout **4.38×** |

B dominates. On one real output file (`hsa-miR-495-3p`, 268,249 rows):
**v1 = 10,289 RefSeq, v2 = 26,219 (2.55×), jaccard 0.214.**

---

## Cause 4 — TargetScan: why it matters

The two parsers answer **different questions**:

- **v1:** *"Which RefSeq transcript **is** the ENST that TargetScan scored?"*
- **v2:** *"Which RefSeq transcripts belong to the **gene** that ENST is part of?"*

TargetScan scores **one transcript's** 3'UTR. Sibling isoforms routinely have entirely
different 3'UTRs — different length, different sites, sometimes no overlap.

> v2 **asserts predicted sites on transcripts that were never evaluated.**
> At transcript granularity this is not noisier — it is wrong.

Of v1's 3,850 missing accessions, **96%** trace to the `9606` filter (Difference A) —
where **v2 is the more correct one**.

---

## Cause 5 — miRmap: a threshold v1 never had ❌

The only tool where v2 predicts **fewer** targets. Controlled test (same 500 UTRs from
v1's own reference, same miRNA): at **engine** level v2 loses nothing — `v1_only = 0`,
and v2 finds **3× more**. The entire shortfall is post-processing:

| # | Cause | Direction | Magnitude |
|---|---|---|---|
| **1** | v2 gates on `ΔG binding ≤ −20`; **v1 has no energy threshold at all** | v2 **much smaller** | **41 → 1 (−97.6%)** |
| **2** | v2 calls `find_targets_with_seed()` bare → miRmap 2 defaults to seeds `[6,7]`; v1 passed `[7,8]` | v2 larger | 41 → 121 |
| **3** | v2 emits every site but the parser reads only the **first** (most-3', not best) | latent | 46% of multi-site transcripts |

**Not the cause:** the 3'UTR reference, the miRmap 1.x→2 migration (costs **zero**
transcripts), or U/T handling.

---

## Cause 5 — miRmap: the ΔG gate is off-scale

The `−20` threshold does not sit at the tail of miRmap's ΔG distribution —
**it sits past the end of it.**

Measured on a full-scale run:

| | value |
|---|--:|
| sites scanned | 242,758 |
| blocks v1 keeps (no threshold) | **84,904** |
| sites passing v2's `ΔG ≤ −20` | **82** |
| ΔG binding — min / median / max | −24.94 / **−7.96** / −3.62 |

The gate discards **~99.9%** of miRmap's output. The `−20` constant appears to have been
copied from miRanda/RNAhybrid — **but it gates a different physical quantity.**

Cause 2 (extra 6mers) partially *masks* Cause 1 — they pull in opposite directions.

---

## Why nobody noticed

The PITA failure is a textbook silent failure — four independent layers of silence:

1. **No error check** — `RNAddG_compute.pl` runs the binary in backticks, ignores `$?`
2. **stderr discarded** — `run_pita()` uses `stderr=devnull`; the segfault message vanishes
3. **Exit status masked** — runs mid-pipeline, so the pipeline reports the *last* command's success
4. **Output looks perfect** — full header, one row per site, `ddG` a plausible number.
   Only 3 unread columns are blank

The others are quieter still: **miRanda, TargetScan and miRmap raise no error at all** —
they just return more, or fewer, rows. The symptom reads as "the tool is permissive"
or "the tool is strict," never as a bug.

---

## What we fixed

| Fix | Where | Status |
|---|---|---|
| Restore miRanda seed filter | `v2/` + `app_v1/parse_result.py` | ✅ Applied |
| Rebuild `RNAddG4` (64-bit, from source) | `tools/PITA64bit/…/Progs/RNAddG4` | ✅ Applied |
| Reproducible build script | `scripts/rebuild_RNAddG4.sh` | ✅ Added |
| PITA accessibility tripwire | `parsePITAResults` | ✅ Applied |
| **TargetScan ENST→RefSeq mapping** | `parse_result.py` | ⏳ **Open** |
| **miRmap ΔG gate / seeds / site pick** | `parse_result.py`, wrapper | ⏳ **Open** |

The tripwire raises if `dGopen` is empty on **all** rows — turning a future silent
recurrence into a loud, self-explaining error.

---

## Verification of the applied fixes

| Check | Result |
|---|---|
| miRanda, miR-21-5p | 358 → **49** = v1 exactly |
| miRanda, miR-1-3p | 144 → **66** = v1 exactly |
| miRanda agreement | **Jaccard 1.000**, zero disagreements |
| Rebuilt `RNAddG4` vs v1 binary | **byte-identical output**, 5/5 sequences |
| PITA end-to-end | `dGopen` 236/236 populated; **target set identical to v1** |
| PITA tripwire | raises on broken output, passes on good |

---

## Decisions we need to make

1. **miRmap ΔG gate — highest priority.** It discards ~99.9% of output. Drop it for v1
   parity, or recalibrate to miRmap's real scale (`≤ −10` sits near the median, keeps ~38%)?
   It must **not** silently reuse miRanda/RNAhybrid's −20.
2. **miRmap seeds** — pass `seed_lengths=[7,8]` explicitly to restore v1? *(Note: 1 and 2
   pull in opposite directions — change together and re-measure.)*
3. **miRmap site selection** — emit/read the best-scoring site rather than the first.
4. **TargetScan mapping** — restore a curated 1:1 ENST→RefSeq table, or keep the symbol
   join? Keeping it is defensible **only** if every consumer collapses to gene symbol anyway.
5. **TargetScan species filter** — v2's `9606` filter is the more correct behaviour. Confirm we keep it.
6. **`RNAddG4` source fix** — an out-of-bounds write (`RNAddG4.c:253-256`) makes it core-dump
   *after* printing. A 2-line bounds guard makes it exit cleanly, predictions unchanged. Apply?
7. **Rebuild the image**, then re-run the full comparison.

---

## Impact on existing results

If any results were generated with v2 as shipped:

| Tool | Status |
|---|---|
| **PITA** | ❌ Invalid — ~20–40× over-predicted, accessibility never applied |
| **miRanda** | ❌ ~7× over-predicted; includes non-seed-matched sites |
| **TargetScan** | ❌ ~2.5× inflated; asserts sites on unevaluated isoforms |
| **miRmap** | ❌ ~99% of real hits discarded by an off-scale threshold |
| **RNAhybrid** | ✅ Trustworthy |

⚠️ **Cross-tool aggregation is the worst affected.** `tool_count`, Venn `sets`/`intersections`
and any consensus logic combine four differently-biased tools — three inflated, one gutted.
Consensus results are not simply noisy; they are systematically skewed.

---

## Appendix — detailed write-ups & repro

**Full root-cause documents** (in repo root):

- `v1_v2_comparison_and_fix.md` — RNAhybrid / miRanda / PITA + the fixes applied
- `v1_v2_comparison_for_targetscan.md` — TargetScan
- `v1_v2_comparison_for_mirmap.md` — miRmap

```bash
# Full comparison harness
cd /bmbl_data/shaohong/projects/playground/v1v2-compare
./run_v1.sh 8 && ./run_v2.sh 8 && python3.11 compare.py

# Is RNAddG4 healthy in an image?
echo "GGGGAAAACCCCAAAAGGGGUUUU" | \
  /opt/PITA64bit/Bin/ViennaRNA/ViennaRNA-1.6/Progs/RNAddG4 -u 0 -s 5 -f 25 -t 25
#   healthy -> prints a tab-separated values line
#   broken  -> prints nothing (exit 139)

# Detect the PITA failure signature in existing output
awk -F'\t' 'NR>1 && $12=="" {c++} END{print "empty dGopen rows:", c+0}' *_pita_results.tab
```
