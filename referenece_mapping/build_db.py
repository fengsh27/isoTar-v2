#!/usr/bin/env python3
"""Build reference_mapping.db from all *_hgnc.tsv files in this directory."""

import os
import sqlite3
import csv

HERE = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(HERE, "reference_mapping.db")

# Map filename prefix -> species code stored in DB
FILE_SPECIES = {
    "cel_WBcel235_3UTRs_hgnc.tsv":      "cel",
    "cfa_CanFam3.1_3UTRs_hgnc.tsv":     "cfa",
    "dme_Release6_3UTRs_hgnc.tsv":      "dme",
    "dre_GRCz11_3UTRs_hgnc.tsv":        "dre",
    "hsa_HG19_only_3UTRs_hgnc.tsv":     "hsa_hg19",
    "hsa_HG38_only_3UTRs_hgnc.tsv":     "hsa_hg38",
    "mdo_MonDom5_3UTRs_hgnc.tsv":       "mdo",
    "mml_Mmul_8.0.1_3UTRs_hgnc.tsv":   "mml",
    "mmu_GRCm38_3UTRs_hgnc.tsv":        "mmu",
    "ptr_Pan_tro3.0_3UTRs_hgnc.tsv":    "ptr",
    "rno_RGSC6_rn6_3UTRs_hgnc.tsv":     "rno",
}

def main():
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
        print("Removed existing {}".format(DB_PATH))

    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()

    cur.execute("""
        CREATE TABLE gene_mapping (
            id      INTEGER PRIMARY KEY AUTOINCREMENT,
            species TEXT NOT NULL,
            raw_id  TEXT NOT NULL,
            symbol  TEXT,
            genename TEXT
        )
    """)

    total = 0
    for filename, species in FILE_SPECIES.items():
        filepath = os.path.join(HERE, filename)
        if not os.path.exists(filepath):
            print("WARNING: {} not found, skipping".format(filepath))
            continue

        rows = []
        with open(filepath, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                rows.append((species, row["raw_ids"], row.get("SYMBOL"), row.get("GENENAME")))

        cur.executemany(
            "INSERT INTO gene_mapping (species, raw_id, symbol, genename) VALUES (?,?,?,?)",
            rows
        )
        print("  {:>10} rows from {} ({})".format(len(rows), filename, species))
        total += len(rows)

    # Indexes for fast lookups
    cur.execute("CREATE INDEX idx_raw_id  ON gene_mapping (raw_id)")
    cur.execute("CREATE INDEX idx_species ON gene_mapping (species)")
    cur.execute("CREATE INDEX idx_symbol  ON gene_mapping (symbol)")

    conn.commit()
    conn.close()
    print("\nDone. {} total rows -> {}".format(total, DB_PATH))

if __name__ == "__main__":
    main()
