#!/usr/bin/env bash
#
# build_targetscan_species_datasets.sh
# ------------------------------------
# Fetch TargetScan's own 3' UTR alignments for the non-human species we support
# and lay them out the way v2/mirna_predicting.py expects, under
# tools/TargetScan/Datasets/<genome>/3utr/targetscan_utr_part_<0..63>.txt.
#
# Dockerfile's `ADD tools /opt/` copies the result into the image. The datasets
# are gitignored (.gitignore: tools/TargetScan/Datasets/*), so this script is
# the reproducible source of truth for them, not the repo.
#
# WHY TargetScan's data and not /opt/reference_files: TargetScan reports its own
# transcript identifiers, and swapping in our reference changes which genes get
# scanned. Measured on hsa-miR-21-5p, the two references agree on only 70% of
# genes (jaccard 0.70) -- because only 33% of shared transcripts even have the
# same 3' UTR length, TargetScan's running ~9% longer. Using their alignments
# keeps TargetScan faithful to what its authors published.
#
# Species covered here are exactly those whose datasets carry an identifier we
# can resolve to RefSeq (see enst_refseq in app_v1/reference_mapping.db):
#     mmu  ENSMUST...   dme  FBtr...   dre  ENSDARG...   cel  <Entrez gene>.<n>
# Human is NOT rebuilt here -- it already lives at Datasets/3utr/ and serves
# both hg19 and hg38.
#
# Worm's UTR id ("171590.0") is its NCBI Entrez gene id plus TargetScan's own
# per-gene counter, so it resolves at gene level, as zebrafish's ENSDARG does.
# TargetScanWorm 6.2 still serves its 5.2 data -- the release notes say worm did
# not change after 5.2 -- which is why its download path says worm_52.
#
# Usage:  scripts/build_targetscan_species_datasets.sh [genome ...]
#         (no args = all four)
# Idempotent: re-running re-downloads and overwrites.
#
# Cost: ~469 MB downloaded, ~4.3 GB on disk / in the image.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST_ROOT="$REPO_ROOT/tools/TargetScan/Datasets"
WORK="${TMPDIR:-/tmp}/ts_species_$$"
NPARTS=64

# genome code -> TargetScan download path
declare -A SRC=(
  [mmu]="mmu_72/mmu_72_data_download/UTR_Sequences.txt.zip"
  [dme]="fly_72/fly_72_data_download/UTR_Sequences.txt.zip"
  [dre]="fish_62/fish_62_data_download/UTR_Sequences.txt.zip"
  [cel]="worm_52/worm_52_data_download/UTR_Sequences.txt.zip"
)

# Columns holding id, species id and sequence, in that order. Worm carries an
# extra "Wormbase Gene ID" column ahead of the species id, so taking 1,4,5 there
# would read the WormBase id as the species and drop the sequence -- a dataset
# that builds cleanly and then silently matches nothing.
declare -A COLS=(
  [mmu]="1,4,5"
  [dme]="1,4,5"
  [dre]="1,4,5"
  [cel]="1,5,6"
)

GENOMES=("$@")
[ ${#GENOMES[@]} -eq 0 ] && GENOMES=(mmu dme dre cel)

mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

for g in "${GENOMES[@]}"; do
  [ -n "${SRC[$g]:-}" ] || { echo "unknown genome: $g" >&2; exit 1; }
  echo "=== $g ==="

  zip="$WORK/$g.zip"
  echo "  downloading..."
  curl -fsS -m 1800 -o "$zip" "https://www.targetscan.org/${SRC[$g]}"

  echo "  extracting..."
  raw="$WORK/${g}_UTR_Sequences.txt"
  unzip -p "$zip" > "$raw"

  # TargetScan ships 5 columns for mouse, fly and zebrafish:
  #   1 transcript/UTR id | 2 gene id | 3 gene symbol | 4 species id | 5 sequence
  # and 6 for worm, with "Wormbase Gene ID" inserted as column 4.
  # targetscan_70.pl wants exactly 3: id, species id, sequence.
  src="$WORK/${g}_3col.txt"
  tail -n +2 "$raw" | cut -f"${COLS[$g]}" > "$src"

  # Guard the column choice: a wrong pick still yields 3 tab-separated fields,
  # so check that the species column is numeric and the sequence is not.
  if ! awk -F'\t' 'NF != 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[-ACGTUNacgtun]+$/ { bad = 1; exit }
                   END { exit bad }' "$src"; then
    echo "  ERROR: columns ${COLS[$g]} did not yield id/species/sequence" >&2
    exit 1
  fi

  # Rows for one gene MUST stay in one part: targetscan_70.pl accumulates rows
  # until the id changes, then scores that set. Splitting a gene across parts
  # would fragment its site groups. The shipped human splitter uses a fixed
  # 37,212 lines, which only works because human is exactly 84 rows per gene;
  # mouse is 52, fly 27, worm 6 and zebrafish 1, so split on the id boundary.
  contiguous=$(cut -f1 "$src" | uniq | wc -l)
  distinct=$(cut -f1 "$src" | sort -u | wc -l)
  if [ "$contiguous" -ne "$distinct" ]; then
    echo "  ERROR: rows are not grouped by id ($contiguous contiguous vs $distinct distinct)" >&2
    exit 1
  fi
  echo "  $(wc -l < "$src") rows, $distinct transcripts"

  dest="$DEST_ROOT/$g/3utr"
  rm -rf "$dest"; mkdir -p "$dest"
  awk -F'\t' -v NP="$NPARTS" -v TG="$distinct" -v OUT="$dest" '
    $1 != prev { gi++; prev = $1 }
    { p = int((gi - 1) * NP / TG); if (p >= NP) p = NP - 1
      print >> (OUT "/targetscan_utr_part_" p ".txt") }
  ' "$src"

  # Every part must exist, even if empty -- mirna_predicting.py loops range(64).
  for i in $(seq 0 $((NPARTS - 1))); do : > /dev/null; [ -f "$dest/targetscan_utr_part_$i.txt" ] || : > "$dest/targetscan_utr_part_$i.txt"; done
  echo "  wrote $(ls "$dest" | wc -l) parts, $(du -sh "$dest" | cut -f1)"
done

echo "done."
