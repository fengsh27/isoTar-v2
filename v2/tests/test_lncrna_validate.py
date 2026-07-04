"""Tests for lncRNA / gene target validation.

Covers the four pieces of the target-validation feature:
  1. normalize_lncrna_id -- version stripping for Ensembl/FlyBase, WormBase intact.
  2. validate_lncrna_targets -- transcript/gene lookup, species scoping.
  3. validate_gene_targets -- symbol + accession lookup (accessions are actually
     checked, unlike resolve_targets).
  4. build_lncrna_db._parse_header -- FASTA header -> (transcript_id, gene).
Plus the POST /api/v1/targets/validate endpoint (needs the Flask app; skipped
where those deps are absent, e.g. outside the backend container).

Run:
    python3 -m unittest v2.tests.test_lncrna_validate
"""

import os
import sys
import sqlite3
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app_v1.lncrna_reference import (  # noqa: E402
    normalize_lncrna_id,
    validate_lncrna_targets,
)
from app_v1.target_resolver import validate_gene_targets  # noqa: E402

# The Flask app pulls flask/celery; guard so pure-logic tests still run locally.
try:
    from app_v1 import app as app_module  # noqa: E402
    _HAVE_APP = True
except Exception:  # pragma: no cover - depends on environment
    app_module = None
    _HAVE_APP = False


def _make_lncrna_db(path, rows):
    """rows: iterable of (species, transcript_id, gene), already normalized."""
    conn = sqlite3.connect(path)
    try:
        c = conn.cursor()
        c.execute("CREATE TABLE lncrna_reference "
                  "(species TEXT NOT NULL, transcript_id TEXT NOT NULL, gene TEXT)")
        c.executemany("INSERT INTO lncrna_reference VALUES (?,?,?)", list(rows))
        conn.commit()
    finally:
        conn.close()


def _make_gene_db(path, rows):
    """rows: iterable of (species, raw_id, symbol, genename)."""
    conn = sqlite3.connect(path)
    try:
        c = conn.cursor()
        c.execute("CREATE TABLE gene_mapping (id INTEGER PRIMARY KEY AUTOINCREMENT, "
                  "species TEXT NOT NULL, raw_id TEXT NOT NULL, symbol TEXT, genename TEXT)")
        c.executemany("INSERT INTO gene_mapping (species, raw_id, symbol, genename) "
                      "VALUES (?,?,?,?)", list(rows))
        conn.commit()
    finally:
        conn.close()


class NormalizeLncrnaIdTests(unittest.TestCase):
    def test_strips_ensembl_transcript_and_gene_versions(self):
        self.assertEqual(normalize_lncrna_id("ENST00000761542.1"), "ENST00000761542")
        self.assertEqual(normalize_lncrna_id("ENSG00000299200.1"), "ENSG00000299200")
        self.assertEqual(normalize_lncrna_id("ENSMUST00000199504.1"), "ENSMUST00000199504")
        self.assertEqual(normalize_lncrna_id("ENSMUSG00000104544.4"), "ENSMUSG00000104544")
        self.assertEqual(normalize_lncrna_id("ENSCAFT00000087990.1"), "ENSCAFT00000087990")

    def test_unversioned_is_unchanged(self):
        self.assertEqual(normalize_lncrna_id("ENST00000761542"), "ENST00000761542")

    def test_flybase_untouched(self):
        self.assertEqual(normalize_lncrna_id("FBtr0347114"), "FBtr0347114")
        self.assertEqual(normalize_lncrna_id("FBgn0267657"), "FBgn0267657")

    def test_wormbase_name_dot_preserved(self):
        # The '.27' is part of the WormBase name, not a version -- must NOT strip.
        self.assertEqual(normalize_lncrna_id("Y51H4A.27"), "Y51H4A.27")
        self.assertEqual(normalize_lncrna_id("WBGene00003296"), "WBGene00003296")

    def test_whitespace_and_empty(self):
        self.assertEqual(normalize_lncrna_id("  ENST00000761542.1  "), "ENST00000761542")
        self.assertIsNone(normalize_lncrna_id(""))
        self.assertIsNone(normalize_lncrna_id(None))
        self.assertIsNone(normalize_lncrna_id("   "))


class ValidateLncrnaTargetsTests(unittest.TestCase):
    def setUp(self):
        fd, self.db = tempfile.mkstemp(suffix=".db", prefix="lnc_ref_")
        os.close(fd)
        # ENSG00000299200 owns two transcripts; a mouse row for species scoping.
        _make_lncrna_db(self.db, [
            ("hsa_hg38", "ENST00000761542", "ENSG00000299200"),
            ("hsa_hg38", "ENST00000761543", "ENSG00000299200"),
            ("mmu", "ENSMUST00000199504", "ENSMUSG00000104544"),
        ])

    def tearDown(self):
        os.remove(self.db)

    def _run(self, targets, genome="hg38"):
        return validate_lncrna_targets(targets, genome, db_path=self.db)

    def test_valid_transcript(self):
        out = self._run(["ENST00000761542"])
        self.assertEqual(out["results"][0]["matched_by"], "transcript")
        self.assertTrue(out["results"][0]["valid"])

    def test_valid_transcript_with_version_suffix(self):
        # User pastes a (different) version; normalization still matches.
        out = self._run(["ENST00000761542.9"])
        self.assertTrue(out["results"][0]["valid"])

    def test_valid_gene_id(self):
        out = self._run(["ENSG00000299200"])
        self.assertEqual(out["results"][0]["matched_by"], "gene")
        self.assertTrue(out["results"][0]["valid"])

    def test_unknown_is_invalid(self):
        out = self._run(["ENST00000000000"])
        self.assertFalse(out["results"][0]["valid"])
        self.assertIsNone(out["results"][0]["matched_by"])

    def test_species_scoping(self):
        # A human transcript must not validate for a mouse job.
        out = self._run(["ENST00000761542"], genome="mmu")
        self.assertFalse(out["results"][0]["valid"])

    def test_batch_order_and_counts(self):
        out = self._run(["ENST00000761542", "bogus", "ENSG00000299200", ""])
        self.assertEqual([r["target"] for r in out["results"]],
                         ["ENST00000761542", "bogus", "ENSG00000299200", ""])
        self.assertEqual(out["valid_count"], 2)
        self.assertEqual(out["invalid"], ["bogus", ""])
        self.assertEqual(out["species"], "hsa_hg38")


class ValidateGeneTargetsTests(unittest.TestCase):
    def setUp(self):
        fd, self.db = tempfile.mkstemp(suffix=".db", prefix="gene_ref_")
        os.close(fd)
        _make_gene_db(self.db, [
            ("hsa_hg38", "NM_000546", "TP53", "tumor protein p53"),
            ("hsa_hg38", "NM_000299", "PKP1", "plakophilin 1"),
            ("mmu", "NM_011640", "Trp53", "transformation related protein 53"),
        ])

    def tearDown(self):
        os.remove(self.db)

    def _run(self, targets, genome="hg38"):
        return validate_gene_targets(targets, genome, ref_db_path=self.db)

    def test_symbol_case_insensitive(self):
        out = self._run(["tp53"])
        self.assertTrue(out["results"][0]["valid"])
        self.assertEqual(out["results"][0]["matched_by"], "symbol")

    def test_accession_with_version(self):
        out = self._run(["NM_000546.7"])
        self.assertTrue(out["results"][0]["valid"])
        self.assertEqual(out["results"][0]["matched_by"], "accession")

    def test_bogus_accession_is_invalid(self):
        # Key difference from resolve_targets: a well-shaped but absent accession
        # is reported invalid, not silently accepted.
        out = self._run(["NM_999999"])
        self.assertFalse(out["results"][0]["valid"])

    def test_unknown_symbol_invalid(self):
        out = self._run(["FOOBAR"])
        self.assertFalse(out["results"][0]["valid"])

    def test_ensembl_ids_are_invalid_for_gene(self):
        # mir-target targets are HGNC symbol / RefSeq only; Ensembl is the
        # lncRNA lane. ENST/ENSG must report invalid on the gene path.
        out = self._run(["ENST00000269305", "ENSG00000141510"])
        self.assertFalse(out["results"][0]["valid"])
        self.assertIsNone(out["results"][0]["matched_by"])
        self.assertFalse(out["results"][1]["valid"])
        self.assertEqual(out["valid_count"], 0)

    def test_species_scoping(self):
        out = self._run(["TP53"], genome="mmu")
        self.assertFalse(out["results"][0]["valid"])


class BuildHeaderParseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, os.path.join(_REPO, "reference_mapping"))
        import build_lncrna_db  # noqa: E402
        cls.parse = staticmethod(build_lncrna_db._parse_header)

    def test_human_header(self):
        h = (">ENST00000761542.1 ncrna scaffold:GRCh38:KI270741.1:44347:130842:1 "
             "gene:ENSG00000299200.1 gene_biotype:lncRNA transcript_biotype:lncRNA")
        self.assertEqual(self.parse(h), ("ENST00000761542", "ENSG00000299200"))

    def test_flybase_header(self):
        h = (">FBtr0347114 ncrna primary_assembly:BDGP6.54:X:1:2:-1 "
             "gene:FBgn0267657 gene_biotype:ncRNA")
        self.assertEqual(self.parse(h), ("FBtr0347114", "FBgn0267657"))

    def test_wormbase_header(self):
        h = ">Y51H4A.27 ncrna chromosome:WBcel235:IV:1:2:1 gene:WBGene00003296"
        self.assertEqual(self.parse(h), ("Y51H4A.27", "WBGene00003296"))

    def test_header_without_gene(self):
        tid, gene = self.parse(">ENST00000761542.1 ncrna description:novel")
        self.assertEqual(tid, "ENST00000761542")
        self.assertIsNone(gene)


@unittest.skipUnless(_HAVE_APP, "flask app deps not available (run in backend container)")
class ValidateEndpointTests(unittest.TestCase):
    """End-to-end against the committed reference DBs via the Flask test client."""

    def setUp(self):
        self.client = app_module.app.test_client()

    def _post(self, body):
        return self.client.post("/api/v1/targets/validate", json=body)

    def test_lncrna_valid_and_invalid(self):
        resp = self._post({
            "target_type": "lncrna", "genome": "hg38",
            "targets": ["ENST00000761542", "ENSG00000299200", "ENST00000000000"],
        })
        self.assertEqual(resp.status_code, 200)
        body = resp.get_json()
        self.assertEqual(body["target_type"], "lncrna")
        self.assertEqual(body["species"], "hsa_hg38")
        by_target = {r["target"]: r for r in body["results"]}
        self.assertTrue(by_target["ENST00000761542"]["valid"])
        self.assertEqual(by_target["ENSG00000299200"]["matched_by"], "gene")
        self.assertFalse(by_target["ENST00000000000"]["valid"])
        self.assertEqual(body["invalid"], ["ENST00000000000"])

    def test_gene_target_id_types(self):
        # mir-target job inputs: gene symbol + RefSeq are accepted; Ensembl is
        # rejected (that's the lncRNA lane). ENST00000269305 maps to TP53 in
        # ensembl_mapping yet must still be invalid here -- we don't silently
        # resolve it, which would let validation disagree with job submission.
        resp = self._post({
            "target_type": "gene", "genome": "hg38",
            "targets": ["TP53", "NM_000546", "ENST00000269305",
                        "ENSG00000141510", "NM_999999"],
        })
        self.assertEqual(resp.status_code, 200)
        by = {r["target"]: r for r in resp.get_json()["results"]}
        self.assertEqual(by["TP53"]["matched_by"], "symbol")
        self.assertEqual(by["NM_000546"]["matched_by"], "accession")
        self.assertFalse(by["ENST00000269305"]["valid"])
        self.assertFalse(by["ENSG00000141510"]["valid"])
        self.assertFalse(by["NM_999999"]["valid"])

    def test_bad_target_type(self):
        resp = self._post({"target_type": "protein", "targets": ["x"]})
        self.assertEqual(resp.status_code, 400)

    def test_targets_must_be_nonempty_list(self):
        self.assertEqual(self._post({"target_type": "lncrna", "targets": []}).status_code, 400)
        self.assertEqual(self._post({"target_type": "lncrna", "targets": "ENST1"}).status_code, 400)


if __name__ == "__main__":
    unittest.main()
