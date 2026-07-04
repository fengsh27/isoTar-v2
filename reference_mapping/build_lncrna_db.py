#!/usr/bin/env python3
"""Build app_v1/lncrna_reference.db from the lncRNA reference FASTAs.

Source of truth is v2/opt/reference_files/*_lncRNAs.fasta (the lncRNA-only
sequence sets the predictor scans). Each FASTA header carries a transcript id
(first token) and a 'gene:<id>' field, e.g.

    >ENST00000761542.1 ncrna ... gene:ENSG00000299200.1 gene_biotype:lncRNA ...

We index (species, transcript_id, gene) with ids version-normalized via the same
normalize_lncrna_id() used at validation time, so unversioned user input matches.
Run manually and commit the resulting .db (as with reference_mapping.db):

    python3 reference_mapping/build_lncrna_db.py
"""

import os
import sys
import sqlite3

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
FASTA_DIR = os.path.join(REPO, "v2", "opt", "reference_files")
DB_PATH = os.path.join(REPO, "app_v1", "lncrna_reference.db")

# Reuse the canonical normalizer so build and validation never diverge.
sys.path.insert(0, REPO)
from app_v1.lncrna_reference import normalize_lncrna_id  # noqa: E402

# Map *_lncRNAs.fasta basename -> species code (same codes as build_db.py).
FILE_SPECIES = {
    "cel_WBcel235_lncRNAs.fasta":    "cel",
    "cfa_CanFam3.1_lncRNAs.fasta":   "cfa",
    "dme_Release6_lncRNAs.fasta":    "dme",
    "dre_GRCz11_lncRNAs.fasta":      "dre",
    "hsa_HG19_lncRNAs.fasta":        "hsa_hg19",
    "hsa_HG38_lncRNAs.fasta":        "hsa_hg38",
    "mdo_MonDom5_lncRNAs.fasta":     "mdo",
    "mml_Mmul_8.0.1_lncRNAs.fasta":  "mml",
    "mmu_GRCm38_lncRNAs.fasta":      "mmu",
    "ptr_Pan_tro3.0_lncRNAs.fasta":  "ptr",
    "rno_RGSC6_rn6_lncRNAs.fasta":   "rno",
}


def _parse_header(line):
    """Return (transcript_id, gene) from a FASTA '>' header, both normalized.

    transcript id is the first whitespace token after '>'; gene is the token
    after 'gene:'. Either may be None if absent."""
    body = line[1:].strip()
    if not body:
        return None, None
    parts = body.split()
    transcript_id = normalize_lncrna_id(parts[0])
    gene = None
    for tok in parts[1:]:
        if tok.startswith("gene:"):
            gene = normalize_lncrna_id(tok[len("gene:"):])
            break
    return transcript_id, gene


def main():
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
        print("Removed existing {}".format(DB_PATH))

    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("""
        CREATE TABLE lncrna_reference (
            species       TEXT NOT NULL,
            transcript_id TEXT NOT NULL,
            gene          TEXT
        )
    """)

    total = 0
    for filename, species in sorted(FILE_SPECIES.items()):
        filepath = os.path.join(FASTA_DIR, filename)
        if not os.path.exists(filepath):
            print("WARNING: {} not found, skipping".format(filepath))
            continue
        rows = []
        with open(filepath, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                transcript_id, gene = _parse_header(line)
                if transcript_id:
                    rows.append((species, transcript_id, gene))
        cur.executemany(
            "INSERT INTO lncrna_reference (species, transcript_id, gene) VALUES (?,?,?)",
            rows,
        )
        print("  {:>8} transcripts from {} ({})".format(len(rows), filename, species))
        total += len(rows)

    cur.execute("CREATE INDEX idx_lnc_tx   ON lncrna_reference (species, transcript_id)")
    cur.execute("CREATE INDEX idx_lnc_gene ON lncrna_reference (species, gene)")
    conn.commit()
    conn.close()
    print("\nDone. lncrna_reference={} rows -> {}".format(total, DB_PATH))


if __name__ == "__main__":
    main()
