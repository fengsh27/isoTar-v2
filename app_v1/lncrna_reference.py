"""Validate user-supplied lncRNA targets against the lncRNA reference index.

The mir-lncrna workflow targets Ensembl/FlyBase/WormBase lncRNA transcripts,
which do not live in reference_mapping.db (that DB is RefSeq/gene-symbol for the
mRNA flow). Valid lncRNA transcript and gene ids come from the reference FASTAs
(v2/opt/reference_files/*_lncRNAs.fasta) and are indexed into lncrna_reference.db
by reference_mapping/build_lncrna_db.py -- schema (species, transcript_id, gene).

Kept Flask-free (like target_resolver.py) so it can be unit-tested and reused by
the build script without importing the app.
"""

import os
import re
import sqlite3

from app_v1.target_resolver import genome_to_species

DEFAULT_LNCRNA_DB = os.environ.get(
    "ISOTAR_LNCRNA_REFERENCE_DB",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "lncrna_reference.db"),
)

# lncRNA target ids across species: Ensembl transcript/gene (ENST/ENSG,
# ENSMUST/ENSMUSG, ENSCAFT/ENSCAFG, ...) and FlyBase (FBtr/FBgn) carry a numeric
# '.<version>' suffix that users usually omit -- strip it so 'ENST00000761542'
# matches the stored 'ENST00000761542.1'. WormBase names (Y51H4A.27, WBGene...)
# have no version and are left intact (the '.27' is part of the name). This is
# the single normalization applied identically at DB-build and validation time.
_LNCRNA_ID_RE = re.compile(r'^(ENS[A-Z]{0,6}[TG]\d+|FB(?:tr|gn)\d+)(?:\.\d+)?')

# SQLite caps a statement at 999 bound parameters; keep IN() batches under that.
_SQL_CHUNK = 900


def normalize_lncrna_id(token):
    """Return a lncRNA transcript OR gene id with any Ensembl/FlyBase version
    suffix removed. WormBase names are returned unchanged. None for empty input."""
    if not token:
        return None
    t = token.strip()
    if not t:
        return None
    m = _LNCRNA_ID_RE.match(t)
    if m:
        return m.group(1)
    return t.split()[0]


def _chunks(seq, size):
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def validate_lncrna_targets(targets, genome, db_path=None):
    """Validate lncRNA target tokens against lncrna_reference for a species.

    Returns a dict:
        {species, results:[{target, valid, matched_by}], valid_count, invalid}
    where matched_by is 'transcript', 'gene' or None. A token is valid when its
    version-normalized form is present as a transcript_id or gene for the species
    (Ensembl/FlyBase namespaces are disjoint, so there is no ambiguity)."""
    db = db_path or DEFAULT_LNCRNA_DB
    species = genome_to_species(genome)

    # Preserve original tokens/order for the response; query each distinct
    # normalized id once.
    order = [(orig, (orig or "").strip()) for orig in targets]
    to_query = []
    seen = set()
    for _orig, tok in order:
        if not tok:
            continue
        nz = normalize_lncrna_id(tok)
        if nz and nz not in seen:
            seen.add(nz)
            to_query.append(nz)

    matched = {}  # normalized id -> 'transcript' | 'gene'
    if to_query and os.path.exists(db):
        conn = sqlite3.connect(db)
        try:
            c = conn.cursor()
            for chunk in _chunks(to_query, _SQL_CHUNK):
                qmarks = ",".join(["?"] * len(chunk))
                c.execute(
                    "SELECT DISTINCT transcript_id FROM lncrna_reference "
                    "WHERE species = ? AND transcript_id IN ({})".format(qmarks),
                    [species] + chunk,
                )
                for (tx,) in c.fetchall():
                    matched[tx] = "transcript"
            for chunk in _chunks(to_query, _SQL_CHUNK):
                qmarks = ",".join(["?"] * len(chunk))
                c.execute(
                    "SELECT DISTINCT gene FROM lncrna_reference "
                    "WHERE species = ? AND gene IN ({})".format(qmarks),
                    [species] + chunk,
                )
                for (g,) in c.fetchall():
                    if g not in matched:  # a transcript hit takes precedence
                        matched[g] = "gene"
        finally:
            conn.close()

    results = []
    invalid = []
    valid_count = 0
    for orig, tok in order:
        mb = matched.get(normalize_lncrna_id(tok)) if tok else None
        valid = mb is not None
        if valid:
            valid_count += 1
        else:
            invalid.append(orig)
        results.append({"target": orig, "valid": valid, "matched_by": mb})

    return {
        "species": species,
        "results": results,
        "valid_count": valid_count,
        "invalid": invalid,
    }
