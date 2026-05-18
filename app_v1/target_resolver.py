"""Resolve user-supplied target lists to RefSeq accessions.

Extracted from app.py so the logic can be tested in isolation (no Flask import).
"""

import os
import re
import sqlite3

DEFAULT_REFERENCE_DB = os.environ.get(
    "ISOTAR_REFERENCE_MAPPING_DB",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference_mapping.db"),
)

# RefSeq accession shape: NM_xxx, NR_xxx, XM_xxx (optional .version)
ACCESSION_RE = re.compile(r"^[A-Z]{2,3}_\d+(\.\d+)?$")


def genome_to_species(genome):
    """Map a UCSC-style genome code to the species value used in gene_mapping."""
    if genome in ("hg19", "hg38"):
        return "hsa_" + genome
    return genome


def resolve_targets(targets, genome, ref_db_path=None):
    """Resolve a mixed list of gene symbols / RefSeq accessions into RefSeq IDs.

    Returns (resolved, unresolved):
        resolved   -- set of unversioned RefSeq accessions (NM_xxx)
        unresolved -- list of input tokens that are not accession-shaped and
                      have no matching row in gene_mapping for the species.
    Accession-shaped tokens are accepted without DB validation.
    """
    db = ref_db_path or DEFAULT_REFERENCE_DB

    resolved = set()
    # Map upper-cased input symbol -> original token; "lose" entries as we find them.
    pending_symbols = {}
    for t in targets:
        t = (t or "").strip()
        if not t:
            continue
        if ACCESSION_RE.match(t):
            resolved.add(t.split(".")[0])
        else:
            pending_symbols.setdefault(t.upper(), t)

    if pending_symbols and os.path.exists(db):
        conn = sqlite3.connect(db)
        try:
            c = conn.cursor()
            qmarks = ",".join(["?"] * len(pending_symbols))
            c.execute(
                "SELECT symbol, raw_id FROM gene_mapping WHERE species = ? "
                "AND symbol COLLATE NOCASE IN ({})".format(qmarks),
                [genome_to_species(genome)] + list(pending_symbols.values()),
            )
            for sym, raw_id in c.fetchall():
                if raw_id:
                    resolved.add(raw_id.split(".")[0])
                pending_symbols.pop(sym.upper(), None)
        finally:
            conn.close()

    unresolved = list(pending_symbols.values())
    return resolved, unresolved
