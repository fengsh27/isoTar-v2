"""Tests for v2.mirna_predicting.filter_utr_fasta.

Run with:
    python3 -m unittest v2.tests.test_filter_utr_fasta
or:
    cd v2/tests && python3 test_filter_utr_fasta.py
"""

import os
import sys
import shutil
import tempfile
import unittest

# Make `v2/` importable when running this file directly.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from v2.mirna_predicting import filter_utr_fasta  # noqa: E402


# Realistic UCSC-style ncbiRefSeq UTR records (accession in header, no gene symbol).
UCSC_FASTA = (
    ">ce11_ncbiRefSeq_NM_059873.7 range=chrI:8377407-8378304 5'pad=0 3'pad=0 strand=- repeatMasking=none\n"
    "ACGTACGTAC\n"
    "ACGTACGTAC\n"
    ">ce11_ncbiRefSeq_NM_182066.7 range=chrI:8377599-8378304 5'pad=0 3'pad=0 strand=- repeatMasking=none\n"
    "TTTTTTTTTT\n"
    ">ce11_ncbiRefSeq_NM_001129046.3 range=chrI:8377601-8378304 5'pad=0 3'pad=0 strand=- repeatMasking=none\n"
    "GGGGGGGGGG\n"
)

# Legacy convention: ">PFX1_PFX2_GENESYM_..." (gene symbol in the 3rd '_' token).
LEGACY_FASTA = (
    ">hsa_HG38_TP53_NM_000546 range=chr17:7565097-7590856 strand=-\n"
    "AAAAAAAAAA\n"
    ">hsa_HG38_BRCA1_NM_007294 range=chr17:41196312-41277500 strand=-\n"
    "CCCCCCCCCC\n"
    ">hsa_HG38_EGFR_NM_005228 range=chr7:55086714-55279321 strand=+\n"
    "GGGGGGGGGG\n"
)


def _headers(path):
    with open(path) as f:
        return [ln.strip() for ln in f if ln.startswith(">")]


def _records(path):
    """Return list of (header, joined-sequence) for sanity-checking sequence preservation."""
    out = []
    header = None
    seq_parts = []
    with open(path) as f:
        for ln in f:
            if ln.startswith(">"):
                if header is not None:
                    out.append((header, "".join(seq_parts)))
                header = ln.strip()
                seq_parts = []
            else:
                seq_parts.append(ln.strip())
    if header is not None:
        out.append((header, "".join(seq_parts)))
    return out


class FilterUtrFastaTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="filter_utr_test_")
        self.ucsc_fa = os.path.join(self.tmp, "ucsc.fa")
        self.legacy_fa = os.path.join(self.tmp, "legacy.fa")
        with open(self.ucsc_fa, "w") as f:
            f.write(UCSC_FASTA)
        with open(self.legacy_fa, "w") as f:
            f.write(LEGACY_FASTA)
        self.out = os.path.join(self.tmp, "out.fa")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    # ---- UCSC headers: match by RefSeq accession ----

    def test_match_by_unversioned_accession(self):
        filter_utr_fasta(self.ucsc_fa, ["NM_059873"], self.out)
        self.assertEqual(len(_headers(self.out)), 1)
        self.assertIn("NM_059873.7", _headers(self.out)[0])

    def test_match_by_versioned_accession_strips_version(self):
        filter_utr_fasta(self.ucsc_fa, ["NM_182066.7"], self.out)
        self.assertEqual(len(_headers(self.out)), 1)
        self.assertIn("NM_182066.7", _headers(self.out)[0])

    def test_versioned_target_matches_differently_versioned_header(self):
        # Target says v9, header is v7 -- version-stripped equality still matches.
        filter_utr_fasta(self.ucsc_fa, ["NM_059873.9"], self.out)
        self.assertEqual(len(_headers(self.out)), 1)

    def test_multiple_targets_keeps_all_matching(self):
        filter_utr_fasta(
            self.ucsc_fa, ["NM_059873", "NM_001129046"], self.out
        )
        headers = _headers(self.out)
        self.assertEqual(len(headers), 2)
        self.assertTrue(any("NM_059873" in h for h in headers))
        self.assertTrue(any("NM_001129046" in h for h in headers))

    def test_no_match_yields_empty_file(self):
        filter_utr_fasta(self.ucsc_fa, ["NM_999999"], self.out)
        self.assertTrue(os.path.exists(self.out))
        self.assertEqual(os.path.getsize(self.out), 0)

    def test_unknown_targets_do_not_error(self):
        filter_utr_fasta(self.ucsc_fa, ["NOT_A_GENE", ""], self.out)
        self.assertEqual(_headers(self.out), [])

    def test_blank_and_whitespace_targets_are_ignored(self):
        filter_utr_fasta(self.ucsc_fa, ["  ", "", "\t", "NM_059873"], self.out)
        self.assertEqual(len(_headers(self.out)), 1)

    # ---- Legacy headers: match by 3rd '_'-token (gene symbol) ----

    def test_legacy_header_matches_by_symbol(self):
        filter_utr_fasta(self.legacy_fa, ["TP53"], self.out)
        headers = _headers(self.out)
        self.assertEqual(len(headers), 1)
        self.assertIn("TP53", headers[0])

    def test_legacy_header_matches_by_embedded_accession(self):
        # Even legacy-style headers contain an NM_ token; accession-matching should also work.
        filter_utr_fasta(self.legacy_fa, ["NM_007294"], self.out)
        headers = _headers(self.out)
        self.assertEqual(len(headers), 1)
        self.assertIn("BRCA1", headers[0])

    def test_mixed_symbol_and_accession_targets(self):
        filter_utr_fasta(
            self.legacy_fa, ["TP53", "NM_005228"], self.out
        )
        headers = _headers(self.out)
        self.assertEqual(len(headers), 2)
        self.assertTrue(any("TP53" in h for h in headers))
        self.assertTrue(any("EGFR" in h for h in headers))

    # ---- Sequence-body preservation ----

    def test_kept_records_preserve_multiline_sequence(self):
        filter_utr_fasta(self.ucsc_fa, ["NM_059873"], self.out)
        recs = _records(self.out)
        self.assertEqual(len(recs), 1)
        # NM_059873 had two 10-character lines -> 20-char concatenated sequence.
        self.assertEqual(recs[0][1], "ACGTACGTAC" * 2)

    def test_unselected_records_excluded_completely(self):
        filter_utr_fasta(self.ucsc_fa, ["NM_059873"], self.out)
        with open(self.out) as f:
            body = f.read()
        self.assertNotIn("TTTTTTTTTT", body)  # NM_182066's sequence
        self.assertNotIn("GGGGGGGGGG", body)  # NM_001129046's sequence

    # ---- Return value ----

    def test_returns_output_path(self):
        ret = filter_utr_fasta(self.ucsc_fa, ["NM_059873"], self.out)
        self.assertEqual(ret, self.out)


if __name__ == "__main__":
    unittest.main()
