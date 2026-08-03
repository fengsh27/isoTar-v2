# isoTar v1 vs v2 — TargetScan Divergence: Root-Cause Analysis

**Date:** 2026-07-21
**Images compared:** `frankfeng78/isotar-v1:1.2.1` vs `frankfeng78/isotar-v2:0.3.18`
**Branch at time of work:** `main` @ `847334e` (chore(release): v0.3.18)
**Companion doc:** [`v1_v2_comparison_and_fix.md`](v1_v2_comparison_and_fix.md) (RNAhybrid / miRanda / PITA),
where TargetScan was listed as "out of scope — different mechanism". This document
retracts that characterisation: the mechanism is *not* different.

---

## TL;DR

TargetScan is the one tool where the 3'UTR reference is **completely irrelevant** — neither
engine feeds it `3utr.fa` or `hsa_HG19_only_3UTRs.fasta`. Both run the **same perl on the
same 6.8 GB pre-aligned dataset with the same miRNA input file**, and therefore produce
**byte-identical raw output**.

**100% of the v1↔v2 divergence happens after the tool runs, inside `parse_result.py`.**
Two independent post-processing differences, pulling in opposite directions:

| # | Difference | Direction | Magnitude |
|---|---|---|---|
| **A** | v2 filters rows to `species_ID == 9606`; v1 does not | v2 **smaller** | 14,124 → 8,426 ENST hits (−40%) |
| **B** | v1 maps ENST→RefSeq **1:1** from a curated table; v2 joins **on gene symbol**, 1:many | v2 **larger** | mean fanout **4.38×** on the hit set |

B dominates. Net on a real single-miRNA output file: **v1 = 10,289 RefSeq, v2 = 26,219
(2.55×), jaccard = 0.214** — reproducing the aggregate figures in the companion doc
(2.0×, 48% of v1 reproduced, jaccard_fair 0.287).

The underlying issue is semantic, not numerical: **v1 reports the transcript TargetScan
actually scored; v2 reports every transcript of that transcript's gene.**

---

## 1. Everything upstream of parsing is identical

TargetScan does **not** consume the per-engine 3'UTR FASTA. Both engines pass it
TargetScan's own pre-split, multi-species aligned UTR set,
`/opt/TargetScan/Datasets/3utr/targetscan_utr_part_{0..63}.txt`, and its matching branch-length
bins. The user/reference UTR file is ignored entirely — `v2/mirna_predicting.py:1128-1130`
says so explicitly, and v1 does the same via `prepareTargetscanTargets()`
(`target.py:87-101`).

Verified identical between the two images:

| item | v1 | v2 |
|---|---|---|
| `TargetScan_70/targetscan_70.pl` | `e756bdb008e456ae593446b0c9823761` | **same** |
| `TargetScan7_BL_PCT/targetscan_70_BL_PCT.pl` | `568bea646adfbbd98664e79d39eb2501` | **same** |
| `Datasets/3utr` | 64 parts, 6.8 GB, `part_0` = `0e9a4e9c1591ffa70f0b9b9498b9a370` | **same** |
| `Datasets/miR_Family_Info.txt` | `65fa07b8cc7d4afcdc2046ceca296f39` | **same** |
| miRNA input file generator | `mirna.py:474-498` | `mirna_predicting.py:468-495` (line-for-line port) |
| invocation | `targetscan_70.pl <mir> <utr> <out>` | identical |

### 1.1 The miRNA input file is provably the same

v1 builds the seed→species map from `miR_Family_Info.txt` (`mirna.py :: getMirFamilyInfo`,
seeding every entry with `'9606'` then appending each species). v2 reads a pre-baked
`miR_Family_Info.json` instead. Rebuilding v1's map from the `.txt` inside the v2 container
and diffing against the `.json`:

```
v1-derived seeds: 4386   json seeds: 4386
identical (order-sensitive): 4386
differing: 0
```

Both then write `<id-without-hsa->\t<seed = seq[1:8], T→U, uppercased>\t<';'-joined species>`.
Same file, same bytes.

### 1.2 Which output file gets parsed differs, but does not matter

- **v1** parses `output_part_%d.txt` — the output of `targetscan_70_BL_PCT.pl`
  (`itutils.py:428-436`, `output_splitted_filename2`; dispatched at `itutils.py:893-897`).
- **v2** parses `<header>_Targetscan_results1.txt` — the merged output of `targetscan_70.pl`
  (`v2/parse_result.py:386`).

BL_PCT is a pass-through annotator: it emits one output row per input row, replacing
columns 12-14 (`Species_in_this_group`, `…_with_this_site_type`, `ORF_overlap`) with
`Branch length score`, `Pct`, `Conserved`. Both files are 14 columns and, on the sample run,
both have 268,249 data rows. Neither parser reads columns 12-14, and the three columns they
*do* read are bit-identical:

```
$ cut -f1,3,9 results1.txt | tail -n +2 | md5sum   ->  f9696790041e84b8e321da7e775455ec
$ cut -f1,3,9 results2.txt | tail -n +2 | md5sum   ->  f9696790041e84b8e321da7e775455ec
```

Applying v1's row filter to either file yields the same 14,124 ENSTs. **Immaterial** — but
worth noting, because it means v2 silently discards the PCT/conservation columns it paid to
compute, and a future conservation filter would have to switch to `results2`.

---

## 2. Difference A — the species filter

TargetScan emits one row per **(transcript, species)** across the whole ortholog alignment.
Critically, column 1 (`a_Gene_ID`) carries the **human** ENST accession on *every* row of the
group, including the mouse and dog rows:

```
a_Gene_ID            species_ID   Site_type
ENST00000000233.5    10090        6mer        <- mouse row, human ENST in col 1
ENST00000000233.5    9615         6mer        <- dog row,   human ENST in col 1
```

Species distribution in the sample file: 47,805 rows for 9606, 46,897 for 9598 (chimp),
44,896 for 9544 (macaque), 35,693 for 9615 (dog), 32,852 for 9913 (cow), …

- **v1** (`itutils.py:834-841`): `len(line) == 14 and line[8] != '6mer'`. **No species test.**
  A human ENST is accepted as a hit when the only non-6mer site in the alignment was found in
  a *non-human* ortholog.
- **v2** (`v2/parse_result.py:145`, `app_v1/parse_result.py:174`):
  `line[2] == "9606" and line[8] != '6mer'`.

| filter | distinct ENST |
|---|--:|
| none | 19,199 |
| v1 (`!= 6mer`) | **14,124** |
| v2 (`9606` + `!= 6mer`) | **8,426** |

v2's ENST set is a **strict subset** of v1's (overlap = 8,426 = all of v2's). This filter is
therefore responsible for essentially all of the "v1 calls missing from v2" bucket — see §4.

---

## 3. Difference B — the ENST→RefSeq mapping (the dominant cause)

Both engines must translate TargetScan's Ensembl transcript IDs into the RefSeq namespace the
rest of the platform uses. They do it in fundamentally different ways.

### 3.1 v1 — a curated 1:1 table

`pred_analysis.py:103` points at `/app/iso/static/public/resources/enst_to_refseq.txt`, a
two-column Biomart dump loaded by `target.py:74-84`:

```python
for line in handler:
    if line[0].strip() != '' and line[1].strip() != '':
        targets_map[line[0].strip()] = line[1].strip()     # plain dict -> last row wins
```

```
rows (non-blank): 42,170
distinct ENST:    38,601
distinct RefSeq:  38,344
ENSTs with >1 RefSeq in the file: 2,884   <- silently collapsed to the last row
```

So: **1 ENST → exactly 1 RefSeq**, and an ENST absent from the table is dropped
(`itutils.py:838`).

### 3.2 v2 — a gene-symbol join

There is no ENST↔RefSeq table in v2. `build_enst_to_refseq_map()`
(`v2/parse_result.py:37-73`, mirrored at `app_v1/parse_result.py:64`) reconstructs one from
`reference_mapping.db` by joining two unrelated tables **on gene symbol**:

```
ensembl_mapping.ensembl_id  ->  symbol  ->  gene_mapping.raw_id  (every RefSeq with that symbol)
```

```
distinct ENST:              130,262
total ENST->RefSeq edges:   745,864
mean fanout:                5.73     max: 368
ENSTs with fanout > 1:      101,632  (78.0%)
distinct RefSeq reachable:  66,814
```

`parseTargetScanResults` then expands every hit through that map
(`v2/parse_result.py:151-160`):

```python
for enst in raw_hits:
    for refseq in enst_to_refseq.get(enst, ()):     # 1 -> N
        ...
```

The optional `targets` argument would clamp the result to RefSeqs present in the run's
`targets.txt`, but that file is absent in the comparison harness (`load_targets_file` returns
`None` → **no filter**), so the full fanout survives.

### 3.3 Isolating its effect

On the actual hit set (7,525 of v2's 8,426 ENSTs are in the map, 89%):

```
mean fanout over mapped hit ENSTs : 4.38
RefSeq before dedup               : 32,996
RefSeq after  dedup               : 26,219
```

Counterfactual — feed **v2's own ENST hits** through **v1's 1:1 map**:

| ENST hit set | through v1's 1:1 map | through v2's symbol map |
|---|--:|--:|
| v1's 14,124 ENSTs | **10,289** | 40,613 |
| v2's 8,426 ENSTs | **6,291** | **26,219** |

Holding the ENST set fixed, the mapping alone inflates results **4.2×** (6,291 → 26,219).
It is the dominant driver, and it more than cancels the 40% reduction from Difference A.

---

## 4. Net effect, measured

Single real output file, `hsa-miR-495-3p` (268,249 data rows), each engine's exact parser
logic applied:

| metric | v1 logic | v2 logic |
|---|--:|--:|
| ENST hits after row filter | 14,124 | 8,426 |
| …of which present in the map | 10,289 (73%) | 7,525 (89%) |
| **final distinct RefSeq** | **10,289** | **26,219** |
| ratio | — | **2.55×** |
| shared | 6,439 | |
| v1-only | 3,850 | |
| v2-only | 19,780 | |
| jaccard | 0.214 | |

Decomposition of the 3,850 **v1-only** accessions:

| origin | count | share |
|---|--:|--:|
| reachable only via ENSTs the `9606` filter dropped (Difference A) | 3,686 | **96%** |
| ENSTs v2 kept, but whose v1 RefSeq is not in v2's symbol fanout | 164 | 4% |

And the 19,780 **v2-only** accessions are, by construction, almost entirely sibling isoforms
introduced by the symbol join.

These figures line up with the companion doc's 5-miRNA aggregate (v1 36,482 / v2 73,371 =
2.0×, 48% of v1 reproduced, jaccard_fair 0.287).

---

## 5. What this actually means

The two parsers answer different questions:

- **v1:** *"Which RefSeq transcript **is** the ENST that TargetScan scored?"*
- **v2:** *"Which RefSeq transcripts belong to the **gene** that ENST is part of?"*

TargetScan scores one specific transcript's 3'UTR. Sibling isoforms of the same gene routinely
have entirely different 3'UTRs — different length, different site content, sometimes no
overlap at all. v2's symbol expansion therefore **asserts predicted sites on transcripts that
were never evaluated**. At transcript granularity that is not merely noisier; it is wrong.

It also skews cross-tool aggregation: TargetScan contributes ~4-6× as many transcript IDs as
the other five tools for the same underlying evidence, which distorts `tool_count`, the Venn
`sets`/`intersections`, and any consensus/intersection logic built on top of them.

The expansion becomes defensible only if every downstream consumer collapses to gene symbol
anyway — in which case the round-trip through RefSeq is pure noise and should be removed
rather than kept.

---

## 6. Recommendations

1. **Replace the symbol join with a real transcript-level mapping.**
   The symbol join is a workaround for a missing table, not a design. Options:
   - lift v1's `enst_to_refseq.txt` out of the v1 image (42,170 rows) — quick, but ~2019-era
     and lossy (2,884 ENSTs silently collapsed to one RefSeq);
   - **preferred:** generate a fresh ENST↔RefSeq table (Biomart, or NCBI/Ensembl
     cross-references) and store it as its own table in `reference_mapping.db`, keeping the
     genuine 1:many relationships that exist at the *transcript* level instead of inventing
     them at the gene level.

2. **Decide the reporting unit deliberately.**
   If the product reports genes (`-tt gene`; the REST result is keyed on `gene_label`), map
   `ENST → symbol` and stop — do not explode to RefSeq and re-collapse. If it reports
   transcripts, item 1 is mandatory.

3. **Settle the species filter on its merits, not by accident.**
   v2's `9606` restriction and v1's conservation-permissive behaviour are both defensible
   positions, but right now the difference is an unreviewed side effect of a rewrite. Whichever
   is chosen should be stated in a comment next to `line[2] == "9606"`.

4. **Parse `results2`, or stop producing it.**
   v2 runs `targetscan_70_BL_PCT.pl` on all 64 parts and merges the output, then parses
   `results1` and never reads it. Either use the `Pct` / `Conserved` columns (a natural place
   for a conservation filter, matching TargetScan's own published behaviour) or drop the
   second pass and reclaim the runtime.

5. **Re-run the v1↔v2 comparison** after items 1-3 to confirm TargetScan converges, alongside
   the PITA/miRanda re-run tracked in the companion doc.

---

## 7. Reproduction

```bash
# --- 1. Confirm the inputs are identical across images ---
for img in frankfeng78/isotar-v1:1.2.1 frankfeng78/isotar-v2:0.3.18; do
  echo "=== $img ==="
  docker run --rm --entrypoint bash "$img" -c '
    md5sum /opt/TargetScan/TargetScan_70/targetscan_70.pl \
           /opt/TargetScan/TargetScan7_BL_PCT/targetscan_70_BL_PCT.pl \
           /opt/TargetScan/Datasets/miR_Family_Info.txt
    ls /opt/TargetScan/Datasets/3utr | wc -l
    md5sum < /opt/TargetScan/Datasets/3utr/targetscan_utr_part_0.txt
    du -sh /opt/TargetScan/Datasets/3utr'
done

# --- 2. Confirm miR_Family_Info.json == v1's .txt-derived map (run inside the v2 image) ---
docker run --rm --entrypoint bash frankfeng78/isotar-v2:0.3.18 -c 'python3.6 - <<PY
import csv, json
m = {}
with open("/opt/TargetScan/Datasets/miR_Family_Info.txt") as f:
    h = csv.reader(f, delimiter="\t"); next(h, None)
    for line in h:
        if len(line) == 7:
            seq, sp = line[1].strip(), line[2].strip()
            if seq and sp:
                m.setdefault(seq, ["9606"])
                if sp not in m[seq]: m[seq].append(sp)
j = json.load(open("/opt/TargetScan/Datasets/miR_Family_Info.json"))
print(len(m), len(j), sum(1 for k in m if list(j.get(k, [])) == m[k]))
PY'

# --- 3. Extract v1'"'"'s mapping table ---
cid=$(docker create frankfeng78/isotar-v1:1.2.1)
docker cp $cid:/app/iso/static/public/resources/enst_to_refseq.txt ./enst_to_refseq.v1.txt
docker rm $cid

# --- 4. Row-filter comparison on a real TargetScan output ---
f=work/mirna_predict_targetscan/Targetscan/hsa-miR-495-3p,WT_Targetscan_results1.txt
awk -F'\t' 'NR>1 && $3=="9606" && $9!="6mer"{g[$1]=1} END{print "v2 filter:", length(g)}' "$f"
awk -F'\t' 'NR>1 &&                $9!="6mer"{g[$1]=1} END{print "v1 filter:", length(g)}' "$f"
```

### 7.1 Full end-to-end attribution script

Run from the repo root, with `enst_to_refseq.v1.txt` from step 3 alongside it:

```python
import csv, re, sqlite3

F = "work/mirna_predict_targetscan/Targetscan/hsa-miR-495-3p,WT_Targetscan_results1.txt"

# v1 map: 1:1, last row wins
v1map = {}
with open("enst_to_refseq.v1.txt") as f:
    h = csv.reader(f, delimiter="\t"); next(h, None)
    for l in h:
        if len(l) > 1 and l[0].strip() and l[1].strip():
            v1map[l[0].strip()] = l[1].strip()

# v2 map: symbol join, 1:many  (mirrors v2/parse_result.py :: build_enst_to_refseq_map)
c = sqlite3.connect("app_v1/reference_mapping.db").cursor()
s2r = {}
for sym, raw in c.execute("SELECT symbol, raw_id FROM gene_mapping "
                          "WHERE species LIKE 'hsa_%' AND symbol IS NOT NULL AND raw_id IS NOT NULL"):
    s2r.setdefault(sym.upper(), set()).add(raw.split(".")[0])
v2map = {}
for eid, sym in c.execute("SELECT ensembl_id, symbol FROM ensembl_mapping "
                          "WHERE ensembl_id IS NOT NULL AND symbol IS NOT NULL"):
    rs = s2r.get(sym.upper())
    if rs:
        v2map.setdefault(eid.split(".")[0], set()).update(rs)

v1e, v2e = set(), set()
with open(F) as f:
    h = csv.reader(f, delimiter="\t"); next(h)
    for line in h:
        if len(line) != 14 or line[8] == "6mer":
            continue
        e = re.sub(r"\.[0-9]+$", "", line[0])
        v1e.add(e)                                   # v1: no species filter
        if line[2] == "9606":
            v2e.add(e)                               # v2: 9606 only

v1r = {v1map[e] for e in v1e if e in v1map}
v2r = set()
for e in v2e:
    v2r |= v2map.get(e, set())

print("ENST hits    v1=%d  v2=%d" % (len(v1e), len(v2e)))
print("final RefSeq v1=%d  v2=%d  (%.2fx)" % (len(v1r), len(v2r), len(v2r) / float(len(v1r))))
print("shared=%d v1_only=%d v2_only=%d jaccard=%.3f" % (
    len(v1r & v2r), len(v1r - v2r), len(v2r - v1r),
    len(v1r & v2r) / float(len(v1r | v2r))))
print("counterfactual: v2's ENSTs through v1's 1:1 map -> %d"
      % len({v1map[e] for e in v2e if e in v1map}))
```

Expected output:

```
ENST hits    v1=14124  v2=8426
final RefSeq v1=10289  v2=26219  (2.55x)
shared=6439 v1_only=3850 v2_only=19780 jaccard=0.214
counterfactual: v2's ENSTs through v1's 1:1 map -> 6291
```

---

## Appendix — key source locations

| item | location |
|---|---|
| v1 TargetScan tool table (cmd / output filenames) | `/app/iso/isotarpkg/itutils.py:428-436` |
| v1 TargetScan runner | `itutils.py :: TargetScanWorker` (269-345) |
| v1 parser (**no species filter**, 1:1 map) | `itutils.py :: parseTargetscanResults` (826-848) |
| v1 mapping table path | `pred_analysis.py:103` → `/app/iso/static/public/resources/enst_to_refseq.txt` |
| v1 mapping table loader (1:1, last-wins) | `target.py :: targetsEnstToRefseqMap` (74-84) |
| v1 miRNA input file writer | `mirna.py:474-498` |
| v1 seed→species map | `mirna.py :: getMirFamilyInfo` (221-245) |
| v2 TargetScan runner | `v2/mirna_predicting.py :: run_targetscan` (411) |
| v2 miRNA input file writer | `v2/mirna_predicting.py :: targetscan_prep` (468-495) |
| v2 TargetScan dataset paths | `v2/mirna_predicting.py:872-873, 1144-1145` |
| v2 **symbol-join** map builder | `v2/parse_result.py :: build_enst_to_refseq_map` (37-73) |
| v2 parser (**9606 filter**, 1:many expansion) | `v2/parse_result.py :: parseTargetScanResults` (131-168) |
| v2 parsed output file (`results1`, not `results2`) | `v2/parse_result.py:386` |
| REST-API copies of both | `app_v1/parse_result.py:64` and `:159` |
