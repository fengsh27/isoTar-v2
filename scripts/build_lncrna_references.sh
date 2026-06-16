#!/usr/bin/env bash
#
# build_lncrna_references.sh
# --------------------------
# Turn the raw Ensembl ncrna dumps in v2/opt/reference_files/lncRNA/*.ncrna.fa.gz
# into filtered, decompressed lncRNA FASTAs that mirror the 3' UTR reference set
# (v2/opt/reference_files/<code>_<assembly>_3UTRs.fasta).
#
# The Ensembl /ncrna/ dumps contain EVERY non-coding biotype (rRNA, snRNA,
# snoRNA, miRNA, piRNA, tRNA, ...). We keep only long-non-coding biotypes. The
# biotype vocabulary differs by release, so the keep-set is a union of every
# spelling Ensembl has used for lncRNA across releases 75-115:
#
#   lncRNA, lincRNA, antisense, antisense_RNA, sense_intronic, sense_overlapping,
#   3prime_overlapping_ncrna, bidirectional_promoter_lncRNA, macro_lncRNA, ncRNA
#
# Notes on the two judgement calls baked in here:
#   * ncRNA (generic): the ONLY home for lncRNAs in the fly (BDGP6.54) and worm
#     (WBcel235) dumps, which lack a proper lncRNA/lincRNA biotype. Included so
#     those species are not empty.
#   * processed_transcript: EXCLUDED (grey-area, often coding-gene-associated;
#     not bona fide lncRNA).
#
# Output names swap "_3UTRs.fasta" for "_lncRNAs.fasta" on the matching 3' UTR
# stem, so v2/mirna_predicting.py can resolve them with the same species-code map.
#
# Usage:  scripts/build_lncrna_references.sh
# Idempotent: re-running overwrites the *_lncRNAs.fasta outputs.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC_DIR="$REPO_ROOT/v2/opt/reference_files/lncRNA"
OUT_DIR="$REPO_ROOT/v2/opt/reference_files"

# Anchored, exact-token match right after "transcript_biotype:" and bounded by a
# space/tab or end-of-line, so e.g. "ncRNA" does not match inside "lncRNA" and
# "antisense" does not swallow "antisense_RNA".
KEEP_RE='transcript_biotype:(lncRNA|lincRNA|antisense|antisense_RNA|sense_intronic|sense_overlapping|3prime_overlapping_ncrna|bidirectional_promoter_lncRNA|macro_lncRNA|ncRNA)([ \t]|$)'

# raw Ensembl filename  ->  output basename (mirrors the 3' UTR stem)
declare -A MAP=(
  [Homo_sapiens.GRCh37.75.ncrna.fa.gz]=hsa_HG19_lncRNAs.fasta
  [Homo_sapiens.GRCh38.ncrna.fa.gz]=hsa_HG38_lncRNAs.fasta
  [Caenorhabditis_elegans.WBcel235.ncrna.fa.gz]=cel_WBcel235_lncRNAs.fasta
  [Canis_lupus_familiaris.CanFam3.1.ncrna.fa.gz]=cfa_CanFam3.1_lncRNAs.fasta
  [Drosophila_melanogaster.BDGP6.54.ncrna.fa.gz]=dme_Release6_lncRNAs.fasta
  [Danio_rerio.GRCz11.ncrna.fa.gz]=dre_GRCz11_lncRNAs.fasta
  [Monodelphis_domestica.monDom5.ncrna.fa.gz]=mdo_MonDom5_lncRNAs.fasta
  [Macaca_mulatta.Mmul_8.0.1.ncrna.fa.gz]=mml_Mmul_8.0.1_lncRNAs.fasta
  [Mus_musculus.GRCm38.ncrna.fa.gz]=mmu_GRCm38_lncRNAs.fasta
  [Pan_troglodytes.Pan_tro_3.0.ncrna.fa.gz]=ptr_Pan_tro3.0_lncRNAs.fasta
  [Rattus_norvegicus.Rnor_6.0.ncrna.fa.gz]=rno_RGSC6_rn6_lncRNAs.fasta
)

printf '%-46s %10s %10s\n' "OUTPUT" "TOTAL" "KEPT(lncRNA)"
printf '%-46s %10s %10s\n' "------" "-----" "------------"

for src in "${!MAP[@]}"; do
  in_path="$SRC_DIR/$src"
  out_path="$OUT_DIR/${MAP[$src]}"
  if [[ ! -f "$in_path" ]]; then
    echo "MISSING source: $in_path" >&2
    exit 1
  fi

  # Keep a header line (and the sequence lines that follow it, until the next
  # header) only when its transcript_biotype is in the keep-set.
  total=$(zcat "$in_path" | grep -c '^>')
  zcat "$in_path" | awk -v re="$KEEP_RE" '
    /^>/ { keep = ($0 ~ re) }
    keep
  ' > "$out_path"
  kept=$(grep -c '^>' "$out_path")
  printf '%-46s %10d %10d\n' "${MAP[$src]}" "$total" "$kept"
done

echo
echo "Done. lncRNA FASTAs written to $OUT_DIR/"
