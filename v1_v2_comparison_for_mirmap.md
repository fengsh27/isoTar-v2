# isoTar v1 vs v2 — miRmap Divergence: Root-Cause Analysis

**Date:** 2026-07-21
**Images compared:** `frankfeng78/isotar-v1:1.2.1` vs `frankfeng78/isotar-v2:0.3.18`
**Branch at time of work:** `main` @ `847334e` (chore(release): v0.3.18)
**Companion docs:** [`v1_v2_comparison_and_fix.md`](v1_v2_comparison_and_fix.md) (RNAhybrid /
miRanda / PITA), [`v1_v2_comparison_for_targetscan.md`](v1_v2_comparison_for_targetscan.md).
Both listed miRmap as "genuinely different software (miRmap 1.x vs miRmap 2) — out of scope".
**This document retracts that.** The software change is real but is *not* the cause of the
divergence.

---

## TL;DR

miRmap is the one tool where v2 predicts **far fewer** targets than v1 — 977 vs 19,285 across
the 5-miRNA benchmark (0.05×, only 2% of v1 reproduced). A controlled experiment (same 500
UTRs from v1's own reference, same miRNA, both engines) shows that at the **engine** level v2
loses nothing at all: every v1 hit is also a v2 hit (`v1_only = 0`), and v2 finds 3× more.

The entire shortfall is introduced afterwards, by three post-processing changes:

| # | Cause | Direction | Magnitude (controlled test) |
|---|---|---|---|
| **1** | v2's parser gates on `ΔG binding <= -20`; **v1 has no energy threshold at all** | v2 **much smaller** | **41 → 1** (‑97.6%) |
| **2** | v2 calls `find_targets_with_seed()` bare → miRmap 2 default seeds `[6,7]`; v1 passed `[7,8]` | v2 larger | 41 → 121 transcripts |
| **3** | v2 emits every site but the parser reads only the **first**, which is the most-3', not the best | latent | 46% of multi-site transcripts affected |

Cause 1 dominates by an order of magnitude and is the answer to "where did the targets go".
Cause 2 partially masks it. Cause 3 is currently invisible but will matter the moment the
threshold is fixed.

**The 3'UTR reference is not the cause** — the controlled test holds it fixed and the
divergence persists at full strength.

---

## 1. The controlled experiment

500 UTRs (200-4,000 nt) taken from **v1's own** `3utr.fa`, scored with **hsa-miR-21-5p**
(`UAGCUUAUCAGACUGAUGUUGA`) by each engine, invoked exactly as its own pipeline invokes it.

| | v1 | v2 |
|---|---|---|
| package | miRmap 1.x, `/opt/miRmap/src` (Python 2.7) | miRmap 2 (`pip`), Python 3.11 |
| entry point | `mirmap.mm(utr, mirna)` + `find_potential_targets_with_seed(...)` | `mirmap.target.find_targets_with_seed(utr, mirna)` |
| seed lengths | **`allowed_lengths=[7,8]`** (explicit) | **`[6,7]`** (library default — not passed) |
| GU wobbles / mismatches | `{7:0, 8:0}` / `{7:0, 8:0}` | no such parameters; exact match only (**equivalent**) |
| sites kept | `take_best=True` → 1 per transcript | **all** sites |

**Raw engine output:**

| metric | count |
|---|--:|
| v1 transcripts with a report block | **41** |
| v2 transcripts with ≥1 site | **121** |
| shared | **41** |
| **v1_only** | **0** |
| v2_only | 80 |

v2's hit set is a **strict superset** of v1's. miRmap 1.x → 2 costs nothing. (v1's 8mer seeds
are still matched by v2 via their embedded 7mer, so dropping `8` from the seed list loses no
transcripts.)

**After each pipeline's parser:**

| parser | kept |
|---|--:|
| v1 (`itutils.py :: parseMirmapResults`) — block exists → keep | **41** |
| v2 (`v2/parse_result.py :: parseMirmapResults`) — `ΔG binding <= -20` | **1** |
| *(hypothetical: v2 gate applied to the **best** site rather than the first)* | *1* |

41 → 1 = **2.4% of v1 retained**, matching the benchmark's reported "2% of v1 reproduced"
almost exactly.

---

## 2. Cause 1 — an energy gate that v1 never had (dominant)

`v2/parse_result.py:264-270`:

```python
matchObj = re.match(r'^\s*ΔG binding \(kcal/mol\)\s+([-+]?\d+\.?\d*)\s*$', lines[i])
if matchObj:
    dg_binding = float(matchObj.group(1))
    # Only add if ΔG binding < -20
    if dg_binding is not None and dg_binding <= -20:
        ...
```

v1 applies **no energy threshold whatsoever** (`itutils.py:735-756`): it keeps a target when a
report block exists — i.e. when the seed search found a site — verified only by the presence
of a numeric line two rows below the `>` header.

### 2.1 v1's report contains no ΔG lines at all

This is not an oversight in v1's parser; the quantity is simply absent. `mirmapRuns`
(`itutils.py:105-135`) evaluates only four features — `eval_tgs_au`, `eval_tgs_pairing3p`,
`eval_tgs_position`, `eval_prob_binomial` — and never `eval_dg_duplex` / `eval_dg_open` /
`eval_dg_binding`. A real v1 report block:

```
1265                           1296
|                              |
ATATGATGCCACTTTGCAGCTAAAATAAGCTCAGTGATACCT
                        |||||||.
          AGTTGTAGTCAGACTATTCGAT
  AU content                     0.47953
  UTR position                   1288.00000
  3' pairing                     0.00000
  Probability (Binomial)         0.31283
```

No `ΔG duplex`, no `ΔG binding`, no `ΔG open`. v2 introduced **both** sides of this gate: the
ΔG lines (synthesized in `_build_mirmap2_block`, `mirna_predicting.py:529-534`) and the
threshold that reads them. There is no v1 behaviour being preserved here.

### 2.2 −20 sits past the end of the distribution

miRmap's `dg_binding` for a seed-anchored site is a small-magnitude quantity. Measured on the
121 v2 hits from the controlled run (first site per transcript):

| statistic | value |
|---|--:|
| min | **−20.52** |
| median | **−9.53** |
| max | −6.63 |

| cutoff | transcripts passing |
|---|--:|
| `<= -5` | 121 / 121 |
| `<= -10` | 46 / 121 |
| `<= -15` | 10 / 121 |
| **`<= -20`** | **1 / 121** |

Independently confirmed at full scale on a real 60 MB miRmap output
(`work/mirna_predict_mirmap/miRmap/hsa-miR-495-3p,WT_miRmap_results.txt`, 242,758 transcripts
scanned):

| metric | count |
|---|--:|
| transcripts with a report block (v1 criterion) | **84,904** (35.0%) |
| passing `ΔG binding <= -20` (v2 criterion) | **82** (0.10% of blocks) |
| ΔG binding min / median / max | −24.94 / −7.96 / −3.62 |
| passing `<= -25` | **0** |

A cutoff of −20 does not filter the distribution; it truncates it just short of the minimum.
**~1 in 1,000 sites survives.**

### 2.3 Where −20 probably came from

−20 is exactly the energy cutoff used for the two neighbouring tools in the same pipeline:
miRanda `-en -20` and RNAhybrid `-e -20`. The most plausible reading is that the constant was
carried across when the miRmap ΔG lines were added.

*This is inference from the value and its context, not proof — no commit message or comment
states it.* But it is a category error regardless of origin: miRanda's and RNAhybrid's −20
gate a **full duplex hybridisation energy**, whereas miRmap's `dg_binding` is a different
model quantity on a different scale (median ≈ −9.5 here). The two are not comparable, and a
threshold calibrated for one is meaningless for the other.

---

## 3. Cause 2 — seed lengths changed silently

v1 (`itutils.py:123-125`) passes its seed policy explicitly:

```python
mim.find_potential_targets_with_seed(
    allowed_lengths=[7,8], allowed_gu_wobbles={7:0,8:0},
    allowed_mismatches={7:0,8:0}, take_best=True)
```

v2 (`mirna_predicting.py:609-610`) passes none of it:

```python
targets = mirmap.target.find_targets_with_seed(utr_sequences[i], mirna_seq_t)
```

miRmap 2's signature is `find_targets_with_seed(host_seq, mirna_seq, seed_start_on_mirna=1,
seed_lengths=None)`, and its body begins:

```python
if seed_lengths is None:
    seed_lengths = [6, 7]
```

So v2 searches **6mer and 7mer** seeds where v1 searched **7mer and 8mer**. Site composition
in the controlled run:

| seed length | sites |
|---|--:|
| 6 | **93** |
| 7 | 41 |
| 8 | — (not searched) |

v2 therefore **adds the weakest, least specific seed class (6mer)** and drops the explicit
strongest one (8mer). This is what takes 41 transcripts to 121, and it pushes v2's counts *up* —
partially masking how severe Cause 1 is. Fixing the threshold without also fixing this would
leave v2 reporting a large body of 6mer-only hits that v1 deliberately excluded.

Note the GU-wobble / mismatch settings need no translation: v1's `{7:0, 8:0}` means zero
wobbles and zero mismatches, and miRmap 2's seed search is a plain reverse-complement regex —
exact match only. The two are already equivalent.

---

## 4. Cause 3 — the parser reads the wrong site (latent)

v1 requested `take_best=True`, so each transcript yields exactly one site and the parser's
fixed offsets are unambiguous.

v2 returns **every** site, and miRmap 2 sorts them by seed end **descending**
(`targets.sort(key=lambda x: x.seed.end, reverse=True)`) — i.e. the site nearest the 3' end of
the UTR comes first, which has nothing to do with score.
`_build_mirmap2_block` (`mirna_predicting.py:497-554`) concatenates all of them, but
`parseMirmapResults` only ever inspects offset `+8` from the `>` header. Subsequent blocks
carry no `>` line, so the loop never revisits them: **every site after the first is invisible
to the parser.**

Each transcript is therefore judged on its most-3' site rather than its best one. In the
controlled run 13 of 121 transcripts had multiple sites, and **6 of those 13 (46%) have
first ≠ best**.

At the current −20 gate this changes nothing (best-site selection also yields 1), which is why
it has gone unnoticed. It becomes a live correctness bug as soon as the threshold is relaxed.

---

## 5. What is *not* the cause

- **The 3'UTR reference.** It does differ (v1 40,374 vs v2 66,869 base accessions, and 62.5%
  of shared accessions differ in length — see the companion doc §3), and it contributes to the
  uncontrolled benchmark. But §1 holds the reference fixed and the divergence persists
  undiminished, so it is not the driver.
- **miRmap 1.x vs miRmap 2.** The migration is real and the report is synthesized rather than
  native, but it costs **zero** transcripts: `v1_only = 0`, and v2's hit set strictly contains
  v1's.
- **The U/T handling.** v1 converts the miRNA to match the target's alphabet
  (`itutils.py:107-117`); v2 unconditionally does `mirna_seq.replace("U", "T")`
  (`mirna_predicting.py:602`). Both end up matching DNA-alphabet UTRs, so on this reference the
  behaviour is the same. (It would diverge on an RNA-alphabet reference — worth a note, not a
  current defect.)

---

## 6. Recommendations

1. **Decide the ΔG gate deliberately — highest priority.** As written it discards ~99.9% of
   miRmap's output. Options:
   - **drop it** for v1 parity (v1 kept every seed-matched transcript); or
   - **recalibrate** to miRmap's actual scale — `<= -10` sits near the median and keeps ~38%.

   Whatever is chosen, it must not silently reuse miRanda/RNAhybrid's −20: that constant gates a
   different physical quantity.

2. **Pass the seed policy explicitly:** `find_targets_with_seed(utr, mirna, seed_lengths=[7,8])`
   to restore v1's behaviour. Leaving it defaulted is how the 6mer class got in.

3. **Fix site selection.** Either have `_build_mirmap2_block` emit the best-scoring site
   (restoring `take_best=True` semantics), or have `parseMirmapResults` scan every block under a
   header instead of only offset `+8`.

4. **Note that 1 and 2 pull in opposite directions.** Relaxing the threshold *and* re-adding
   6mers would overshoot badly. Change them together and re-measure against v1 before
   concluding anything.

5. **Re-run the v1↔v2 comparison** after 1-3, alongside the PITA/miRanda/TargetScan re-runs
   tracked in the companion docs.

---

## 7. Reproduction

Both scripts below live in the session scratchpad, not in the repo. `$SP` is any scratch
directory.

```bash
# --- Extract v1's reference and take a 500-UTR subset ---
cid=$(docker create frankfeng78/isotar-v1:1.2.1)
docker cp $cid:/app/iso/static/public/resources/3utr.fa "$SP/v1_3utr.fa"
docker rm $cid

python3 - <<EOF
out=[]; n=0; hdr=None; seq=[]
for line in open("$SP/v1_3utr.fa"):
    if line.startswith(">"):
        if hdr and 200 <= len("".join(seq)) <= 4000:
            out.append((hdr,"".join(seq))); n+=1
            if n>=500: break
        hdr=line.strip(); seq=[]
    else: seq.append(line.strip())
with open("$SP/sub500.fa","w") as f:
    for h,s in out: f.write(h+"\n"+s+"\n")
EOF
```

### 7.1 v1 engine (miRmap 1.x, Python 2.7)

> **Gotcha:** the script contains a literal `Δ`, so Python 2.7 requires a PEP 263 header
> (`# -*- coding: utf-8 -*-`) or it dies with
> `SyntaxError: Non-ASCII character '\xce'`. Same rule CLAUDE.md documents for `v2/*.py`.

```python
# -*- coding: utf-8 -*-
import re, mirmap, json
MIR = "UAGCUUAUCAGACUGAUGUUGA"
hdr=None; seq=[]; recs=[]
for line in open("/w/sub500.fa"):
    if line.startswith(">"):
        if hdr: recs.append((hdr,"".join(seq)))
        hdr=line.strip(); seq=[]
    else: seq.append(line.strip())
if hdr: recs.append((hdr,"".join(seq)))

res={}
for h,s in recs:
    acc = re.search(r'(N[MR]_\d+)', h).group(1)
    mirna = MIR
    if s.find('U') == -1 and mirna.find('U') != -1:
        mirna = mirna.replace('U','T')
    try:
        mim = mirmap.mm(s, mirna)
        mim.find_potential_targets_with_seed(allowed_lengths=[7,8],
            allowed_gu_wobbles={7:0,8:0}, allowed_mismatches={7:0,8:0}, take_best=True)
        mim.end_sites
        mim.eval_tgs_au(with_correction=False)
        mim.eval_tgs_pairing3p(with_correction=False)
        mim.eval_tgs_position(with_correction=False)
        mim.eval_prob_binomial()
        rep = mim.report()
    except ValueError:
        rep = ""
    if rep.strip():
        res[acc] = True
json.dump(res, open("/w/v1_mirmap.json","w"))
print("v1: %d/%d transcripts produced a report block" % (len(res), len(recs)))
```

```bash
docker run --rm -v "$SP:/w" --entrypoint bash frankfeng78/isotar-v1:1.2.1 \
  -c 'cd /w && python /w/v1_run.py'
# -> v1: 41/500 transcripts produced a report block
```

### 7.2 v2 engine (miRmap 2, Python 3.11)

```python
import re, json
import mirmap.target, mirmap.thermo
MIR = "UAGCUUAUCAGACUGAUGUUGA".replace("U","T")
hdr=None; seq=[]; recs=[]
for line in open("/w/sub500.fa"):
    if line.startswith(">"):
        if hdr: recs.append((hdr,"".join(seq)))
        hdr=line.strip(); seq=[]
    else: seq.append(line.strip())
if hdr: recs.append((hdr,"".join(seq)))

res={}
for h,s in recs:
    acc = re.search(r'(N[MR]_\d+)', h).group(1)
    targets = mirmap.target.find_targets_with_seed(s, MIR)   # bare, as v2 calls it
    if not targets: continue
    sites=[]
    for t in targets:
        try:   d = mirmap.thermo.calc_dg_duplex(t)["dg_binding"]
        except Exception: d = None
        sites.append({"dg_binding": d, "seed_len": t.seed_length, "seed_end": t.seed.end})
    dgs=[x["dg_binding"] for x in sites if x["dg_binding"] is not None]
    res[acc] = {"n_sites": len(sites), "sites": sites,
                "first_dg": sites[0]["dg_binding"],
                "best_dg": min(dgs) if dgs else None}
json.dump(res, open("/w/v2_mirmap.json","w"))
print("v2: %d/%d transcripts produced >=1 site" % (len(res), len(recs)))
```

```bash
docker run --rm -v "$SP:/w" --entrypoint bash frankfeng78/isotar-v2:0.3.18 \
  -c 'cd /w && python3.11 /w/v2_run.py'
# -> v2: 121/500 transcripts produced >=1 site
```

### 7.3 Confirm miRmap 2's seed default

```bash
docker run --rm --entrypoint bash frankfeng78/isotar-v2:0.3.18 -c \
 'python3.11 -c "import inspect, mirmap.target
print(inspect.signature(mirmap.target.find_targets_with_seed))
print(inspect.getsource(mirmap.target.find_targets_with_seed)[:900])"'
# -> (host_seq, mirna_seq, seed_start_on_mirna=1, seed_lengths=None)
# -> if seed_lengths is None: seed_lengths = [6, 7]
```

### 7.4 Threshold effect on an existing full-scale output

```bash
python3 - <<'EOF'
import io, re, statistics
f = "work/mirna_predict_mirmap/miRmap/hsa-miR-495-3p,WT_miRmap_results.txt"
lines = io.open(f, encoding='utf-8').readlines()
tot = blocks = 0; dgs = []
for i, l in enumerate(lines):
    if not re.match(r'^>[^,]+,.*?\s+(\S+)\s*$', l): continue
    tot += 1
    if i+2 < len(lines) and re.match(r'.*[0-9]+.*', lines[i+2]):
        blocks += 1
        if i+8 < len(lines):
            m = re.match(r'^\s*ΔG binding \(kcal/mol\)\s+([-+]?\d+\.?\d*)', lines[i+8])
            if m: dgs.append(float(m.group(1)))
print("scanned=%d  blocks(v1 keeps)=%d  dG<=-20 (v2 keeps)=%d"
      % (tot, blocks, sum(1 for d in dgs if d <= -20)))
print("dG binding min=%.2f median=%.2f max=%.2f"
      % (min(dgs), statistics.median(dgs), max(dgs)))
EOF
# -> scanned=242758  blocks(v1 keeps)=84904  dG<=-20 (v2 keeps)=82
# -> dG binding min=-24.94 median=-7.96 max=-3.62
```

---

## Appendix — key source locations

| item | location |
|---|---|
| v1 miRmap invocation (seeds `[7,8]`, `take_best`) | `/app/iso/isotarpkg/itutils.py:105-135` (`mirmapRuns`), call at `:123-125` |
| v1 feature evaluation (no ΔG evaluated) | `itutils.py:127-130` |
| v1 per-target report writer | `itutils.py :: miRmapWorker` (366-392), header + report at `:378-380` |
| v1 parser (**no energy threshold**) | `itutils.py :: parseMirmapResults` (735-756) |
| v1 miRmap package | `/opt/miRmap/src/mirmap/__init__.py` (1.x, Python 2.7) |
| v2 miRmap runner | `v2/mirna_predicting.py :: run_mirmap` (557-618) |
| v2 seed search (**bare defaults → `[6,7]`**) | `v2/mirna_predicting.py:609-610` |
| v2 synthesized 1.x-style report block | `v2/mirna_predicting.py :: _build_mirmap2_block` (497-554); ΔG lines at `:529-534` |
| v2 shard failure tolerance | `v2/mirna_predicting.py:443-453` |
| v2 parser (**`ΔG binding <= -20`**) | `v2/parse_result.py :: parseMirmapResults` (245-275); threshold at `:268` |
| v2 parser offsets (`+2` numeric, `+8` ΔG binding) | `v2/parse_result.py:260-264` |
| miRmap 2 seed default `[6,7]` | `mirmap/target.py :: find_targets_with_seed` (site-packages, Python 3.11) |
