"""Tests for app_v1.target_resolver.resolve_targets.

Run with:
    python3 -m unittest v2.tests.test_target_resolver
"""

import os
import sys
import shutil
import sqlite3
import tempfile
import unittest

# Make the repo root importable when running this file directly.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app_v1.target_resolver import (  # noqa: E402
    ACCESSION_RE,
    genome_to_species,
    resolve_targets,
)


def _make_ref_db(path, gene_rows):
    """Build a minimal reference_mapping.db with the supplied gene_mapping rows.

    gene_rows: iterable of (species, raw_id, symbol, genename)
    """
    conn = sqlite3.connect(path)
    try:
        c = conn.cursor()
        c.execute("""
            CREATE TABLE gene_mapping (
                id      INTEGER PRIMARY KEY AUTOINCREMENT,
                species TEXT NOT NULL,
                raw_id  TEXT NOT NULL,
                symbol  TEXT,
                genename TEXT
            )
        """)
        c.executemany(
            "INSERT INTO gene_mapping (species, raw_id, symbol, genename) VALUES (?,?,?,?)",
            gene_rows,
        )
        conn.commit()
    finally:
        conn.close()


class GenomeToSpeciesTests(unittest.TestCase):
    def test_human_genomes_get_hsa_prefix(self):
        self.assertEqual(genome_to_species("hg19"), "hsa_hg19")
        self.assertEqual(genome_to_species("hg38"), "hsa_hg38")

    def test_other_species_passthrough(self):
        for code in ("cel", "cfa", "dme", "dre", "mdo", "mml", "mmu", "ptr", "rno"):
            self.assertEqual(genome_to_species(code), code)


class AccessionRegexTests(unittest.TestCase):
    def test_matches_nm_with_and_without_version(self):
        self.assertTrue(ACCESSION_RE.match("NM_000546"))
        self.assertTrue(ACCESSION_RE.match("NM_000546.6"))

    def test_matches_nr_and_xm(self):
        self.assertTrue(ACCESSION_RE.match("NR_046018"))
        self.assertTrue(ACCESSION_RE.match("XM_011517506"))

    def test_does_not_match_symbols_or_junk(self):
        for s in ("TP53", "tp53", "ENST00000269305", "NM_", "NMXXXX", "123_456", ""):
            self.assertIsNone(ACCESSION_RE.match(s), s)


class ResolveTargetsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="resolve_targets_test_")
        self.db = os.path.join(self.tmp, "ref.db")
        _make_ref_db(
            self.db,
            gene_rows=[
                ("hsa_hg38", "NM_000546",    "TP53",     "tumor protein p53"),
                ("hsa_hg38", "NM_001126112", "TP53",     "tumor protein p53"),
                ("hsa_hg38", "NM_001126113", "TP53",     "tumor protein p53"),
                ("hsa_hg38", "NM_007294",    "BRCA1",    "breast cancer 1"),
                ("hsa_hg38", "NM_005228",    "EGFR",     "EGFR"),
                ("hsa_hg38", "NM_018672",    "C1orf141", "chr1 ORF 141"),
                ("hsa_hg19", "NM_hg19_only", "ONLYHG19", "only hg19"),
                ("mmu",      "Trp53_mouse",  "Trp53",    "mouse trp53"),
            ],
        )

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    # ---- Accession-shaped inputs ----

    def test_accession_passthrough_strips_version(self):
        resolved, unresolved = resolve_targets(["NM_000546.6"], "hg38", ref_db_path=self.db)
        self.assertEqual(resolved, {"NM_000546"})
        self.assertEqual(unresolved, [])

    def test_unversioned_accession_kept_as_is(self):
        resolved, unresolved = resolve_targets(["NM_000546"], "hg38", ref_db_path=self.db)
        self.assertEqual(resolved, {"NM_000546"})
        self.assertEqual(unresolved, [])

    def test_accession_not_validated_against_db(self):
        # Bogus accession passes regex and is accepted unvalidated — spec says it's
        # the user's responsibility for accession typos.
        resolved, unresolved = resolve_targets(["NM_99999999"], "hg38", ref_db_path=self.db)
        self.assertEqual(resolved, {"NM_99999999"})
        self.assertEqual(unresolved, [])

    # ---- Symbol lookups ----

    def test_symbol_resolves_to_all_refseqs(self):
        resolved, unresolved = resolve_targets(["TP53"], "hg38", ref_db_path=self.db)
        self.assertEqual(resolved, {"NM_000546", "NM_001126112", "NM_001126113"})
        self.assertEqual(unresolved, [])

    def test_symbol_lookup_is_case_insensitive(self):
        resolved, unresolved = resolve_targets(["tp53"], "hg38", ref_db_path=self.db)
        self.assertEqual(resolved, {"NM_000546", "NM_001126112", "NM_001126113"})
        self.assertEqual(unresolved, [])

    def test_mixed_case_symbol_in_db_matches_uppercase_input(self):
        # gene_mapping stores "C1orf141" — input "C1ORF141" should still resolve.
        resolved, unresolved = resolve_targets(["C1ORF141"], "hg38", ref_db_path=self.db)
        self.assertEqual(resolved, {"NM_018672"})
        self.assertEqual(unresolved, [])

    def test_unknown_symbol_returns_in_unresolved(self):
        resolved, unresolved = resolve_targets(["T53"], "hg38", ref_db_path=self.db)
        self.assertEqual(resolved, set())
        self.assertEqual(unresolved, ["T53"])

    def test_unresolved_preserves_original_token_case(self):
        resolved, unresolved = resolve_targets(["foobar"], "hg38", ref_db_path=self.db)
        self.assertEqual(unresolved, ["foobar"])

    def test_partial_resolution(self):
        # TP53 + BRCA1 resolve; T53 (typo) does not.
        resolved, unresolved = resolve_targets(
            ["TP53", "T53", "BRCA1"], "hg38", ref_db_path=self.db,
        )
        self.assertEqual(
            resolved,
            {"NM_000546", "NM_001126112", "NM_001126113", "NM_007294"},
        )
        self.assertEqual(unresolved, ["T53"])

    # ---- Mixed inputs ----

    def test_mixed_symbol_and_accession(self):
        resolved, unresolved = resolve_targets(
            ["TP53", "NM_007294", "EGFR"], "hg38", ref_db_path=self.db,
        )
        self.assertEqual(
            resolved,
            {"NM_000546", "NM_001126112", "NM_001126113", "NM_007294", "NM_005228"},
        )
        self.assertEqual(unresolved, [])

    # ---- Species filtering ----

    def test_species_filter_hg38(self):
        # ONLYHG19 has rows only in hsa_hg19; should not resolve for hg38.
        resolved, unresolved = resolve_targets(["ONLYHG19"], "hg38", ref_db_path=self.db)
        self.assertEqual(resolved, set())
        self.assertEqual(unresolved, ["ONLYHG19"])

    def test_species_filter_hg19(self):
        resolved, unresolved = resolve_targets(["ONLYHG19"], "hg19", ref_db_path=self.db)
        self.assertEqual(resolved, {"NM_hg19_only"})
        self.assertEqual(unresolved, [])

    def test_human_symbol_does_not_match_mouse_row(self):
        # TP53 is human-only; mouse uses Trp53. Querying hg38 for Trp53 should miss.
        resolved, unresolved = resolve_targets(["Trp53"], "hg38", ref_db_path=self.db)
        self.assertEqual(resolved, set())
        self.assertEqual(unresolved, ["Trp53"])

    def test_mouse_species_resolves_mouse_symbol(self):
        resolved, unresolved = resolve_targets(["Trp53"], "mmu", ref_db_path=self.db)
        self.assertEqual(resolved, {"Trp53_mouse"})
        self.assertEqual(unresolved, [])

    # ---- Edge cases ----

    def test_empty_list(self):
        resolved, unresolved = resolve_targets([], "hg38", ref_db_path=self.db)
        self.assertEqual(resolved, set())
        self.assertEqual(unresolved, [])

    def test_blank_and_whitespace_tokens_ignored(self):
        resolved, unresolved = resolve_targets(
            ["", "  ", "\t", "TP53"], "hg38", ref_db_path=self.db,
        )
        self.assertEqual(
            resolved,
            {"NM_000546", "NM_001126112", "NM_001126113"},
        )
        self.assertEqual(unresolved, [])

    def test_none_token_is_ignored(self):
        resolved, unresolved = resolve_targets(
            [None, "TP53"], "hg38", ref_db_path=self.db,
        )
        self.assertEqual(
            resolved,
            {"NM_000546", "NM_001126112", "NM_001126113"},
        )
        self.assertEqual(unresolved, [])

    def test_duplicate_symbols_dedup_in_unresolved(self):
        # Same unknown symbol appearing twice should only show once in unresolved.
        resolved, unresolved = resolve_targets(
            ["FOOBAR", "FOOBAR"], "hg38", ref_db_path=self.db,
        )
        self.assertEqual(resolved, set())
        self.assertEqual(unresolved, ["FOOBAR"])

    def test_missing_db_marks_symbols_unresolved(self):
        # Without a DB to consult, only accession-shaped tokens can pass through.
        resolved, unresolved = resolve_targets(
            ["TP53", "NM_000546"], "hg38", ref_db_path="/no/such.db",
        )
        self.assertEqual(resolved, {"NM_000546"})
        self.assertEqual(unresolved, ["TP53"])


if __name__ == "__main__":
    unittest.main()
