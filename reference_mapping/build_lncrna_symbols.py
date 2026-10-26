#!/usr/bin/env python3
"""Build app_v1/lncrna_symbols.tsv (species, gene, symbol) for literature verification.

Source: an Ensembl ncRNA FASTA dump whose headers carry `gene:` and `gene_symbol:`
fields (GRCh38 does; the GRCh37.75 dump has no gene_symbol, so human hg19 reuses
the GRCh38 symbols through the assembly-stable ENSG ids). Only long non-coding
biotypes are kept -- the same keep-set as scripts/build_lncrna_references.sh.

Usage:
    python3 reference_mapping/build_lncrna_symbols.py <species> <ncrna.fa[.gz]> [--out PATH]

    python3 reference_mapping/build_lncrna_symbols.py hsa Homo_sapiens.GRCh38.ncrna.fa.gz

Rows are appended per species: an existing file keeps the other species' rows.
The verification code maps hsa_hg19 / hsa_hg38 onto the species key "hsa".
"""
import argparse
import csv
import gzip
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DEFAULT_OUT = os.path.join(REPO, "app_v1", "lncrna_symbols.tsv")

sys.path.insert(0, REPO)
from app_v1.lncrna_reference import normalize_lncrna_id  # noqa: E402

KEEP = {"lncRNA", "lincRNA", "antisense", "antisense_RNA", "sense_intronic", "sense_overlapping",
        "3prime_overlapping_ncrna", "bidirectional_promoter_lncRNA", "macro_lncRNA", "ncRNA"}
FIELD = re.compile(r"(\S+?):(\S+)")


def parse_header(line):
    fields = {}
    for tok in line[1:].split()[1:]:
        m = FIELD.match(tok)
        if m:
            fields.setdefault(m.group(1), m.group(2))
    return fields


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("species")
    ap.add_argument("fasta")
    ap.add_argument("--out", default=DEFAULT_OUT)
    a = ap.parse_args()

    opener = gzip.open if a.fasta.endswith(".gz") else open
    symbols = {}
    with opener(a.fasta, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            fields = parse_header(line)
            if fields.get("transcript_biotype") not in KEEP:
                continue
            gene, symbol = fields.get("gene"), fields.get("gene_symbol")
            if gene and symbol:
                symbols.setdefault(normalize_lncrna_id(gene), symbol)

    rows = []
    if os.path.exists(a.out):
        with open(a.out, newline="", encoding="utf-8") as f:
            rows = [r for r in csv.DictReader(f, delimiter="\t") if r["species"] != a.species]
    rows += [{"species": a.species, "gene": g, "symbol": s} for g, s in sorted(symbols.items())]
    with open(a.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["species", "gene", "symbol"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("%s: %d lncRNA genes with a symbol -> %s (%d rows total)" % (a.species, len(symbols), a.out, len(rows)))


if __name__ == "__main__":
    main()
