"""Tests for running TargetScan on worm (genome code `cel`).

Worm's TargetScan UTR id ("171687.1") is its NCBI Entrez gene id plus
TargetScan's own per-gene counter. The counter is stripped like a version
suffix, and the Entrez id resolves to RefSeq through the enst_refseq table under
build WBcel235 -- at gene level, as zebrafish's ENSDARG ids already do.

The tests against the shipped reference_mapping.db skip if it is absent.

Run with:
    python3 -m unittest v2.tests.test_targetscan_worm
"""

import json
import os
import shutil
import sqlite3
import sys
import tempfile
import unittest

try:
    from unittest import mock
except ImportError:  # pragma: no cover
    mock = None

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

import v2.mirna_predicting as runner  # noqa: E402
import v2.parse_result as parser_v2  # noqa: E402
import app_v1.parse_result as parser_app  # noqa: E402

_SHIPPED_DB = os.path.join(_REPO, "app_v1", "reference_mapping.db")


def _make_enst_refseq_db(path, rows):
    """rows: iterable of (build, enst, refseq)."""
    conn = sqlite3.connect(path)
    try:
        conn.execute("CREATE TABLE enst_refseq (build TEXT NOT NULL, "
                     "enst TEXT NOT NULL, refseq TEXT NOT NULL)")
        conn.executemany("INSERT INTO enst_refseq VALUES (?,?,?)", rows)
        conn.commit()
    finally:
        conn.close()


def _ts_row(utr_id, species, site="7mer-m8"):
    r = [""] * 14
    r[0] = utr_id
    r[2] = species
    r[8] = site
    return r


class GenomeConfigAgreementTests(unittest.TestCase):
    """The runner decides which genomes TargetScan runs for; the parsers decide
    which genomes the API accepts and how hits are mapped. If they drift, a job
    is either rejected for no reason or runs and maps nothing."""

    def test_cel_is_supported_everywhere(self):
        self.assertTrue(runner.targetscan_supported("cel"))
        self.assertTrue(parser_v2.targetscan_supported("cel"))
        self.assertTrue(parser_app.targetscan_supported("cel"))

    def test_runner_and_parsers_agree_on_genomes(self):
        runner_genomes = set(runner._TARGETSCAN_GENOME_DIR)
        self.assertEqual(runner_genomes, set(parser_v2.TARGETSCAN_GENOMES))
        self.assertEqual(runner_genomes, set(parser_app.TARGETSCAN_GENOMES))

    def test_runner_and_parsers_agree_on_taxids(self):
        for genome in runner._TARGETSCAN_GENOME_DIR:
            self.assertEqual(runner._GENOME_TAXID[genome],
                             parser_v2.genome_taxid(genome), genome)
            self.assertEqual(runner._GENOME_TAXID[genome],
                             parser_app.genome_taxid(genome), genome)

    def test_cel_taxid_and_build(self):
        # 6239 is C. elegans; the worm alignment also carries 6238 (C. briggsae)
        # and four more nematodes, whose rows must not count as worm targets.
        self.assertEqual(parser_v2.genome_taxid("cel"), "6239")
        self.assertEqual(parser_v2._GENOME_BUILD["cel"], "WBcel235")
        self.assertEqual(runner.targetscan_utr_dir("cel"),
                         os.path.join(runner.TARGETSCAN, "Datasets", "cel/3utr"))

    def test_cel_has_no_precomputed_bins(self):
        self.assertIsNone(runner.targetscan_bins_dir("cel"))


class WormMappingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="ts_worm_test_")
        self.db = os.path.join(self.tmp, "ref.db")
        _make_enst_refseq_db(self.db, [
            ("WBcel235", "171687", "NM_001380927"),
            ("WBcel235", "171687", "NM_001380932"),
            ("WBcel235", "171726", "NM_058473"),
            # Same key under another build must not leak into worm.
            ("GRCz11", "171726", "NM_zebrafish_only"),
        ])
        self.out = os.path.join(self.tmp, "ts.txt")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, rows):
        with open(self.out, "w") as f:
            f.write("\t".join(["a"] * 14) + "\n")
            for r in rows:
                f.write("\t".join(r) + "\n")

    def test_map_reads_only_the_worm_build(self):
        mp = parser_v2.build_enst_to_refseq_map(self.db, genome="cel")
        self.assertEqual(mp, {
            "171687": {"NM_001380927", "NM_001380932"},
            "171726": {"NM_058473"},
        })

    def test_worm_hits_resolve_to_refseq(self):
        # 171687.1 and 171687.2 are two UTRs of one gene with the same site; they
        # must collapse, not double-count. 6238 (C. briggsae) and 6mer rows are
        # dropped.
        self._write([
            _ts_row("171687.1", "6239"),
            _ts_row("171687.2", "6239"),
            _ts_row("171726.0", "6239", site="7mer-1a"),
            _ts_row("999999.0", "6238"),
            _ts_row("888888.0", "6239", site="6mer"),
        ])
        for parser in (parser_v2, parser_app):
            mp = parser.build_enst_to_refseq_map(self.db, genome="cel")
            out = parser.parseTargetScanResults(
                self.out, {}, enst_to_refseq=mp,
                species_id=parser.genome_taxid("cel"))
            self.assertEqual(
                sorted(out["prediction"]["TargetScan"]),
                ["NM_001380927", "NM_001380932", "NM_058473"],
                parser.__name__)

    def test_human_filter_finds_nothing_in_worm_output(self):
        # The silent-failure guard: parsing worm output at 9606 returns an empty
        # list rather than an error, which is why the genome must reach here.
        self._write([_ts_row("171687.1", "6239")])
        mp = parser_v2.build_enst_to_refseq_map(self.db, genome="cel")
        out = parser_v2.parseTargetScanResults(self.out, {}, enst_to_refseq=mp,
                                               species_id="9606")
        self.assertEqual(out["prediction"]["TargetScan"], [])


@unittest.skipUnless(mock is not None, "unittest.mock unavailable")
class WormTargetScanPrepTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="ts_worm_prep_")
        self.mirfam = os.path.join(self.tmp, "miR_Family_Info.json")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _species_written(self, family_info):
        with open(self.mirfam, "w") as f:
            json.dump(family_info, f)
        with mock.patch.object(runner, "targetscan_mirfam_path",
                               return_value=self.mirfam):
            runner.targetscan_prep("UGAGGUAGUAGGUUGUAUAGUU",
                                   "cel-let-7-5p", self.tmp, genome="cel")
        path = os.path.join(self.tmp, "cel-let-7-5p_targetscan.txt")
        with open(path) as f:
            return f.readline().rstrip("\n").split("\t")

    def test_worm_taxon_added_when_seed_matches_vertebrate_family(self):
        # let-7's seed is in the shipped vertebrate file; without 6239 appended,
        # targetscan_70.pl would skip every worm UTR and return nothing.
        fields = self._species_written({"GAGGUAG": ["9606", "10090"]})
        self.assertEqual(fields[1], "GAGGUAG")
        self.assertIn("6239", fields[2].split(";"))

    def test_worm_taxon_used_when_seed_unknown(self):
        fields = self._species_written({})
        self.assertEqual(fields[2], "6239")


@unittest.skipUnless(os.path.exists(_SHIPPED_DB), "reference_mapping.db not present")
class ShippedWormMapTests(unittest.TestCase):
    def test_shipped_db_carries_worm_rows(self):
        mp = parser_v2.build_enst_to_refseq_map(_SHIPPED_DB, genome="cel")
        self.assertGreater(len(mp), 19000)
        self.assertIn("NM_001380927", mp["171687"])
        self.assertEqual(mp["171726"], {"NM_058473"})
        for refseqs in mp.values():
            for r in refseqs:
                self.assertTrue(r.startswith("NM_"), r)


if __name__ == "__main__":
    unittest.main()
