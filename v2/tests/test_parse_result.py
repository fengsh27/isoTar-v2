"""Tests for v2.parse_result helpers added for TargetScan post-processing:
    - load_targets_file
    - build_enst_to_refseq_map
    - parseTargetScanResults (ENST->RefSeq conversion + targets filter)

Run with:
    python3 -m unittest v2.tests.test_parse_result
"""

import os
import sys
import shutil
import sqlite3
import tempfile
import unittest

# Make `v2/` importable when running this file directly.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from v2.parse_result import (  # noqa: E402
    build_enst_to_refseq_map,
    load_targets_file,
    parseTargetScanResults,
    parseDMISOResults,
    _mirna_seed_match_pattern,
)

# The live result-building path (app_v1) carries the lncRNA-aware parsers; the
# v2 copy is the legacy gene-only version. Import the app_v1 PITA parser and the
# lenient lncRNA id extractor directly to lock the header-row regression.
from app_v1.parse_result import (  # noqa: E402
    parsePITAResults as parsePITAResults_app,
    _extract_lncrna_transcript_id,
    _extract_transcript_id,
)


def _make_ref_db(path, gene_rows, ensembl_rows):
    """Build a minimal reference_mapping.db with the supplied rows.

    gene_rows:    iterable of (species, raw_id, symbol, genename)
    ensembl_rows: iterable of (ensembl_id, symbol)
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
        c.execute("CREATE TABLE ensembl_mapping (ensembl_id TEXT, symbol TEXT)")
        c.executemany(
            "INSERT INTO gene_mapping (species, raw_id, symbol, genename) VALUES (?,?,?,?)",
            gene_rows,
        )
        c.executemany(
            "INSERT INTO ensembl_mapping (ensembl_id, symbol) VALUES (?,?)",
            ensembl_rows,
        )
        conn.commit()
    finally:
        conn.close()


class LoadTargetsFileTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="targets_test_")
        self.path = os.path.join(self.tmp, "targets.txt")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_missing_file_returns_none(self):
        self.assertIsNone(load_targets_file(self.path))

    def test_none_arg_returns_none(self):
        self.assertIsNone(load_targets_file(None))

    def test_empty_string_arg_returns_none(self):
        self.assertIsNone(load_targets_file(""))

    def test_empty_file_returns_none(self):
        open(self.path, "w").close()
        self.assertIsNone(load_targets_file(self.path))

    def test_keeps_only_refseq_shapes(self):
        with open(self.path, "w") as f:
            f.write("NM_000546\n")
            f.write("TP53\n")          # bare symbol -- dropped
            f.write("FOOBAR\n")        # garbage -- dropped
            f.write("NM_001126112\n")
            f.write("\n")
            f.write("  \n")
        self.assertEqual(
            load_targets_file(self.path),
            {"NM_000546", "NM_001126112"},
        )

    def test_strips_version_suffix(self):
        with open(self.path, "w") as f:
            f.write("NM_000546.6\n")
            f.write("NM_001126112.3\n")
        self.assertEqual(
            load_targets_file(self.path),
            {"NM_000546", "NM_001126112"},
        )

    def test_accepts_xm_and_nr_prefixes(self):
        # Regex is [A-Z]{2,3}_\d+ -- should accept XM_ (predicted mRNA) and NR_ (non-coding RNA).
        with open(self.path, "w") as f:
            f.write("NM_000546\n")
            f.write("XM_011517506\n")
            f.write("NR_046018\n")
        self.assertEqual(
            load_targets_file(self.path),
            {"NM_000546", "XM_011517506", "NR_046018"},
        )


class BuildEnstToRefseqMapTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="enstmap_test_")
        self.db = os.path.join(self.tmp, "ref.db")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_missing_db_returns_empty(self):
        self.assertEqual(build_enst_to_refseq_map("/no/such/path.db"), {})

    def test_basic_one_to_many_join(self):
        _make_ref_db(
            self.db,
            gene_rows=[
                ("hsa_hg38", "NM_000546",    "TP53",  "tumor protein p53"),
                ("hsa_hg38", "NM_001126112", "TP53",  "tumor protein p53"),
                ("hsa_hg38", "NM_007294",    "BRCA1", "breast cancer 1"),
            ],
            ensembl_rows=[
                ("ENST00000269305", "TP53"),
                ("ENST00000357654", "BRCA1"),
            ],
        )
        mp = build_enst_to_refseq_map(self.db)
        self.assertEqual(mp["ENST00000269305"], {"NM_000546", "NM_001126112"})
        self.assertEqual(mp["ENST00000357654"], {"NM_007294"})

    def test_case_insensitive_symbol_match(self):
        # gene_mapping stores "C1orf141" (mixed case) but ensembl_mapping has "C1ORF141".
        _make_ref_db(
            self.db,
            gene_rows=[("hsa_hg38", "NM_018672", "C1orf141", "chr1 ORF 141")],
            ensembl_rows=[("ENST00000425777", "C1ORF141")],
        )
        mp = build_enst_to_refseq_map(self.db)
        self.assertEqual(mp.get("ENST00000425777"), {"NM_018672"})

    def test_species_filter_excludes_nonhuman(self):
        # Same ENST appears with both a human symbol (resolves) and a mouse one (filtered out).
        _make_ref_db(
            self.db,
            gene_rows=[
                ("hsa_hg38", "NM_000546",       "TP53",  "..."),
                ("mmu",      "NR_mouse_id",     "Trp53", "..."),
            ],
            ensembl_rows=[
                ("ENST00000269305", "TP53"),
                ("ENST00000269305", "Trp53"),
            ],
        )
        mp = build_enst_to_refseq_map(self.db)
        self.assertEqual(mp["ENST00000269305"], {"NM_000546"})
        # The mouse raw_id must NOT leak into the map anywhere.
        for refseqs in mp.values():
            self.assertNotIn("NR_mouse_id", refseqs)

    def test_strips_versions_from_both_ids(self):
        _make_ref_db(
            self.db,
            gene_rows=[("hsa_hg38", "NM_000546.6", "TP53", "...")],
            ensembl_rows=[("ENST00000269305.10", "TP53")],
        )
        mp = build_enst_to_refseq_map(self.db)
        self.assertIn("ENST00000269305", mp)
        self.assertNotIn("ENST00000269305.10", mp)
        self.assertEqual(mp["ENST00000269305"], {"NM_000546"})

    def test_orphan_enst_dropped(self):
        # ENST whose symbol has no row in gene_mapping should not appear.
        _make_ref_db(
            self.db,
            gene_rows=[("hsa_hg38", "NM_000546", "TP53", "...")],
            ensembl_rows=[
                ("ENST00000269305", "TP53"),
                ("ENST00000000001", "MYSTERY_GENE"),
            ],
        )
        mp = build_enst_to_refseq_map(self.db)
        self.assertIn("ENST00000269305", mp)
        self.assertNotIn("ENST00000000001", mp)

    def test_includes_both_hg19_and_hg38(self):
        # species LIKE 'hsa_%' -- both hg19 and hg38 rows should join.
        _make_ref_db(
            self.db,
            gene_rows=[
                ("hsa_hg19", "NM_hg19_id", "TP53", "..."),
                ("hsa_hg38", "NM_hg38_id", "TP53", "..."),
            ],
            ensembl_rows=[("ENST00000269305", "TP53")],
        )
        mp = build_enst_to_refseq_map(self.db)
        self.assertEqual(mp["ENST00000269305"], {"NM_hg19_id", "NM_hg38_id"})


class ParseTargetScanResultsTests(unittest.TestCase):
    """TargetScan output format: TSV, 14 columns.
        col 0 = transcript ID
        col 2 = species id ('9606' == human)
        col 8 = site type ('6mer' is excluded)
    parseTargetScanResults skips rows that don't match these constraints."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="parse_ts_test_")
        self.out = os.path.join(self.tmp, "ts.txt")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write_rows(self, rows):
        # Header (14 cols) + each data row (14 cols).
        with open(self.out, "w") as f:
            header = ["a"] * 14
            f.write("\t".join(header) + "\n")
            for row in rows:
                assert len(row) == 14, "test row must have 14 cols"
                f.write("\t".join(row) + "\n")

    def _row(self, tid, species="9606", site="7mer-m8"):
        r = [""] * 14
        r[0] = tid
        r[2] = species
        r[8] = site
        return r

    def test_passthrough_when_no_map(self):
        self._write_rows([
            self._row("ENST00000269305.10"),
            self._row("ENST00000357654.4", site="8mer"),
        ])
        out = parseTargetScanResults(self.out, {})
        self.assertEqual(
            out["prediction"]["TargetScan"],
            ["ENST00000269305", "ENST00000357654"],
        )

    def test_converts_enst_to_refseq_when_map_given(self):
        enst_map = {
            "ENST00000269305": {"NM_000546", "NM_001126112"},
            "ENST00000357654": {"NM_007294"},
        }
        self._write_rows([
            self._row("ENST00000269305"),
            self._row("ENST00000357654", site="8mer"),
        ])
        out = parseTargetScanResults(self.out, {}, enst_to_refseq=enst_map)
        self.assertEqual(
            set(out["prediction"]["TargetScan"]),
            {"NM_000546", "NM_001126112", "NM_007294"},
        )

    def test_filters_by_targets_set(self):
        enst_map = {
            "ENST00000269305": {"NM_000546", "NM_001126112"},
            "ENST00000357654": {"NM_007294"},
        }
        targets = {"NM_000546", "NM_007294"}  # exclude NM_001126112
        self._write_rows([
            self._row("ENST00000269305"),
            self._row("ENST00000357654", site="8mer"),
        ])
        out = parseTargetScanResults(self.out, {}, enst_to_refseq=enst_map, targets=targets)
        self.assertEqual(
            set(out["prediction"]["TargetScan"]),
            {"NM_000546", "NM_007294"},
        )

    def test_empty_targets_set_drops_everything(self):
        # An empty (but non-None) targets set means "filter against nothing".
        enst_map = {"ENST00000269305": {"NM_000546"}}
        self._write_rows([self._row("ENST00000269305")])
        out = parseTargetScanResults(self.out, {}, enst_to_refseq=enst_map, targets=set())
        self.assertEqual(out["prediction"]["TargetScan"], [])

    def test_drops_non_human_rows(self):
        enst_map = {
            "ENST00000269305":    {"NM_000546"},
            "ENSMUST00000010000": {"NM_mouse_x"},
        }
        self._write_rows([
            self._row("ENST00000269305"),
            self._row("ENSMUST00000010000", species="10090"),
        ])
        out = parseTargetScanResults(self.out, {}, enst_to_refseq=enst_map)
        self.assertEqual(out["prediction"]["TargetScan"], ["NM_000546"])

    def test_drops_6mer_sites(self):
        enst_map = {
            "ENST00000269305": {"NM_000546"},
            "ENST00000357654": {"NM_007294"},
        }
        self._write_rows([
            self._row("ENST00000269305", site="7mer-m8"),
            self._row("ENST00000357654", site="6mer"),
        ])
        out = parseTargetScanResults(self.out, {}, enst_to_refseq=enst_map)
        self.assertEqual(out["prediction"]["TargetScan"], ["NM_000546"])

    def test_missing_file_yields_empty_list(self):
        out = parseTargetScanResults("/no/such/file.txt", {})
        self.assertEqual(out["prediction"]["TargetScan"], [])

    def test_drops_enst_not_in_map(self):
        # ENSTs absent from the conversion map are silently dropped.
        enst_map = {"ENST00000269305": {"NM_000546"}}
        self._write_rows([
            self._row("ENST00000269305"),
            self._row("ENST99999999999"),
        ])
        out = parseTargetScanResults(self.out, {}, enst_to_refseq=enst_map)
        self.assertEqual(out["prediction"]["TargetScan"], ["NM_000546"])

    def test_dedup_when_multiple_ensts_share_a_refseq(self):
        enst_map = {
            "ENST_A": {"NM_X"},
            "ENST_B": {"NM_X"},  # same NM appears under two ENSTs
        }
        self._write_rows([
            self._row("ENST_A"),
            self._row("ENST_B"),
        ])
        out = parseTargetScanResults(self.out, {}, enst_to_refseq=enst_map)
        self.assertEqual(out["prediction"]["TargetScan"], ["NM_X"])

    def test_result_dict_is_not_overwritten(self):
        # Other keys must be preserved when we set prediction.TargetScan.
        self._write_rows([self._row("ENST00000269305")])
        rd = {"header": "hsa-let-7a", "type": "WT"}
        out = parseTargetScanResults(self.out, rd)
        self.assertEqual(out["header"], "hsa-let-7a")
        self.assertEqual(out["type"], "WT")
        self.assertIn("TargetScan", out["prediction"])


class SeedPatternTests(unittest.TestCase):
    """_mirna_seed_match_pattern: reverse-complement of miRNA positions 2-8."""

    def test_known_mirna(self):
        # hsa-miR-495-3p: seed (pos 2-8) reverse-complement is GTTTGTT
        self.assertEqual(
            _mirna_seed_match_pattern("AAACAAACAUGGUGCACUUCUU"), "GTTTGTT"
        )

    def test_rna_u_and_lowercase_normalized(self):
        self.assertEqual(
            _mirna_seed_match_pattern("aaacaaacauggugcacuucuu"), "GTTTGTT"
        )

    def test_returns_seven_mer(self):
        self.assertEqual(len(_mirna_seed_match_pattern("AAACAAACAUGGUGCACUUCUU")), 7)

    def test_poly_a(self):
        # rev of AAAAAAAA is AAAAAAAA; seed slice [-8:-1] -> AAAAAAA; complement TTTTTTT
        self.assertEqual(_mirna_seed_match_pattern("AAAAAAAA"), "TTTTTTT")

    def test_exactly_eight_nt_ok(self):
        self.assertEqual(len(_mirna_seed_match_pattern("ACGUACGU")), 7)

    def test_none_returns_empty(self):
        self.assertEqual(_mirna_seed_match_pattern(None), "")

    def test_empty_returns_empty(self):
        self.assertEqual(_mirna_seed_match_pattern(""), "")

    def test_too_short_returns_empty(self):
        # < 8 nt -> '' -> caller applies no filter
        self.assertEqual(_mirna_seed_match_pattern("AAACAAA"), "")


class ParseDMISOResultsTests(unittest.TestCase):
    """parseDMISOResults: keep only targets whose sequence contains the seed."""

    # miR-495-3p seed-match pattern is GTTTGTT (see SeedPatternTests)
    MIRNA = "AAACAAACAUGGUGCACUUCUU"

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="dmiso_test_")
        self.path = os.path.join(self.tmp, "x_DMISO_results.txt")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, rows):
        with open(self.path, "w") as f:
            f.write("Target ID\tTarget Sequence\tPrediction Score\n")
            for r in rows:
                f.write("\t".join(r) + "\n")

    def test_keeps_seed_match_drops_non_match(self):
        self._write([
            ("hg19_refGene_NM_144722 range=chr5:1-9", "TTAGGTTTGTTGGATG", "0.991"),  # has GTTTGTT
            ("hg19_refGene_NM_018034 range=chr5:2-9", "AGAATCTCATTTGAGAG", "0.995"),  # no match
        ])
        rd = parseDMISOResults(self.path, {}, mirna_sequence=self.MIRNA)
        self.assertEqual(rd["prediction"]["DMISO"], ["NM_144722"])

    def test_no_filter_when_sequence_missing(self):
        self._write([
            ("hg19_refGene_NM_144722 range=x", "AAAAAAAAAA", "0.991"),
            ("hg19_refGene_NM_018034 range=x", "CCCCCCCCCC", "0.995"),
        ])
        rd = parseDMISOResults(self.path, {}, mirna_sequence=None)
        self.assertEqual(
            sorted(rd["prediction"]["DMISO"]), ["NM_018034", "NM_144722"]
        )

    def test_dedup(self):
        self._write([
            ("hg19_refGene_NM_144722 range=a", "GTTTGTTAAA", "0.991"),
            ("hg19_refGene_NM_144722 range=b", "CCGTTTGTTC", "0.992"),
        ])
        rd = parseDMISOResults(self.path, {}, mirna_sequence=self.MIRNA)
        self.assertEqual(rd["prediction"]["DMISO"], ["NM_144722"])

    def test_case_insensitive_target_sequence(self):
        self._write([
            ("hg19_refGene_NM_144722 range=a", "ttagctgtttgttgg", "0.991"),  # lowercase match
        ])
        rd = parseDMISOResults(self.path, {}, mirna_sequence=self.MIRNA)
        self.assertEqual(rd["prediction"]["DMISO"], ["NM_144722"])

    def test_missing_file_yields_empty_list(self):
        rd = parseDMISOResults(os.path.join(self.tmp, "nope.txt"), {},
                               mirna_sequence=self.MIRNA)
        self.assertEqual(rd["prediction"]["DMISO"], [])

    def test_short_rows_skipped(self):
        with open(self.path, "w") as f:
            f.write("Target ID\tTarget Sequence\tPrediction Score\n")
            f.write("only_one_column\n")
            f.write("hg19_refGene_NM_144722 range=a\tGTTTGTTAAA\t0.991\n")
        rd = parseDMISOResults(self.path, {}, mirna_sequence=self.MIRNA)
        self.assertEqual(rd["prediction"]["DMISO"], ["NM_144722"])

    def test_preserves_existing_result_dict_keys(self):
        self._write([("hg19_refGene_NM_144722 range=a", "GTTTGTTAAA", "0.991")])
        rd = parseDMISOResults(
            self.path, {"header": "x", "type": "WT"}, mirna_sequence=self.MIRNA
        )
        self.assertEqual(rd["header"], "x")
        self.assertEqual(rd["type"], "WT")
        self.assertIn("DMISO", rd["prediction"])


class ParsePITAResultsHeaderTests(unittest.TestCase):
    """Regression for the lncRNA PITA crash: the PITA output starts with a
    header row (``UTR ... ddG``). The gene flow skips it because
    _extract_transcript_id returns None for the "UTR" cell, but the lenient
    lncRNA extractor returns "UTR" verbatim, which used to reach float("ddG")
    and raise ValueError -- crashing the whole result build."""

    HEADER = "UTR\tmicroRNA\tStart\tEnd\tSeed\tLoop\tdGduplex\tdG5\tdG3\tdG0\tdG1\tdGopen\tddG\n"

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="pita_test_")
        self.path = os.path.join(self.tmp, "x_PITA_results.tab")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _row(self, target, ddg):
        # 13 columns; only col0 (target) and col12 (ddG) are read.
        return "{}\thsa-let-7a-3p,WT\t1\t8\t7:0:0\t0\t-5\t-5\t0\t\t\t\t{}\n".format(target, ddg)

    def _write(self, rows):
        with open(self.path, "w") as f:
            f.write(self.HEADER)
            for target, ddg in rows:
                f.write(self._row(target, ddg))

    def test_lncrna_header_row_does_not_crash_and_is_excluded(self):
        hdr = "ENST00000357591.3 ncrna chromosome:GRCh38:4:1:2:1 gene_biotype:lncRNA"
        self._write([(hdr, "-11.9"), ("ENST00000111111.1 ncrna foo", "-6.62")])
        rd = parsePITAResults_app(self.path, {}, id_extractor=_extract_lncrna_transcript_id)
        # ddG <= -10 kept (full header collapsed to the version-stripped
        # accession, matching every other lncRNA parser); ddG > -10 dropped.
        self.assertEqual(rd["prediction"]["PITA"], ["ENST00000357591"])
        self.assertNotIn("UTR", rd["prediction"]["PITA"])

    def test_gene_flow_unchanged(self):
        self._write([("hg19_refGene_NM_000051 range=a", "-11.9"),
                     ("hg19_refGene_NM_000546 range=b", "-6.62")])
        rd = parsePITAResults_app(self.path, {})  # default gene extractor
        self.assertEqual(rd["prediction"]["PITA"], ["NM_000051"])

    def test_non_numeric_ddg_row_skipped(self):
        self._write([("ENST00000357591.3 x", "n/a"), ("ENST00000222222.2 y", "-12.0")])
        rd = parsePITAResults_app(self.path, {}, id_extractor=_extract_lncrna_transcript_id)
        self.assertEqual(rd["prediction"]["PITA"], ["ENST00000222222"])


if __name__ == "__main__":
    unittest.main()
