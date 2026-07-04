"""Tests for v2.mirna_predicting.filter_lncrna_fasta.

The lncRNA target filter keeps a reference record if its transcript id OR its
'gene:' id (version-normalized) is in the target set -- so a gene-id target fans
out to all that gene's transcripts. Normalization must match the backend
validator (normalize_lncrna_id).

Run with:
    python3 -m unittest v2.tests.test_filter_lncrna_fasta
"""

import os
import sys
import shutil
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from v2.mirna_predicting import filter_lncrna_fasta  # noqa: E402


# Two transcripts share gene ENSG00000299200; a third gene; plus FlyBase and
# WormBase records to exercise the non-Ensembl id shapes.
LNCRNA_FASTA = (
    ">ENST00000761542.1 ncrna scaffold:GRCh38 gene:ENSG00000299200.1 "
    "gene_biotype:lncRNA transcript_biotype:lncRNA\n"
    "ACGTACGTAC\n"
    "ACGTACGTAC\n"
    ">ENST00000761543.1 ncrna scaffold:GRCh38 gene:ENSG00000299200.1 "
    "gene_biotype:lncRNA transcript_biotype:lncRNA\n"
    "TTTTTTTTTT\n"
    ">ENST00000999999.2 ncrna scaffold:GRCh38 gene:ENSG00000123456.4 "
    "gene_biotype:lncRNA transcript_biotype:lncRNA\n"
    "GGGGGGGGGG\n"
    ">FBtr0347114 ncrna primary_assembly:BDGP6.54 gene:FBgn0267657 "
    "gene_biotype:ncRNA\n"
    "CCCCCCCCCC\n"
    ">Y51H4A.27 ncrna chromosome:WBcel235 gene:WBGene00003296 gene_biotype:ncRNA\n"
    "AAAAAAAAAA\n"
)


def _headers(path):
    with open(path) as f:
        return [ln.strip() for ln in f if ln.startswith(">")]


class FilterLncrnaFastaTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="filter_lncrna_test_")
        self.fa = os.path.join(self.tmp, "lncrna.fa")
        with open(self.fa, "w") as f:
            f.write(LNCRNA_FASTA)
        self.out = os.path.join(self.tmp, "out.fa")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, targets):
        filter_lncrna_fasta(self.fa, targets, self.out)
        return _headers(self.out)

    def test_match_single_transcript(self):
        headers = self._run(["ENST00000761542"])
        self.assertEqual(len(headers), 1)
        self.assertIn("ENST00000761542.1", headers[0])

    def test_transcript_version_insensitive(self):
        # User pastes a different version; version-stripped equality still matches.
        headers = self._run(["ENST00000761542.9"])
        self.assertEqual(len(headers), 1)

    def test_gene_id_fans_out_to_all_transcripts(self):
        headers = self._run(["ENSG00000299200"])
        self.assertEqual(len(headers), 2)
        self.assertTrue(all("gene:ENSG00000299200" in h for h in headers))
        self.assertTrue(any("ENST00000761542" in h for h in headers))
        self.assertTrue(any("ENST00000761543" in h for h in headers))

    def test_gene_id_with_version(self):
        headers = self._run(["ENSG00000299200.7"])
        self.assertEqual(len(headers), 2)

    def test_flybase_transcript_and_gene(self):
        self.assertEqual(len(self._run(["FBtr0347114"])), 1)
        self.assertEqual(len(self._run(["FBgn0267657"])), 1)

    def test_wormbase_name_preserved(self):
        # The '.27' is part of the WormBase name, not a version.
        headers = self._run(["Y51H4A.27"])
        self.assertEqual(len(headers), 1)
        self.assertIn("Y51H4A.27", headers[0])

    def test_wormbase_gene_id(self):
        self.assertEqual(len(self._run(["WBGene00003296"])), 1)

    def test_unknown_target_keeps_nothing(self):
        filter_lncrna_fasta(self.fa, ["ENST00000000000"], self.out)
        self.assertEqual(os.path.getsize(self.out), 0)

    def test_blank_targets_ignored(self):
        headers = self._run(["", "  ", "\t", "ENST00000761542"])
        self.assertEqual(len(headers), 1)

    def test_mixed_transcript_and_gene_targets(self):
        # One transcript from gene A + gene B's id -> that transcript + B's transcript.
        headers = self._run(["ENST00000761542", "ENSG00000123456"])
        self.assertEqual(len(headers), 2)
        self.assertTrue(any("ENST00000761542" in h for h in headers))
        self.assertTrue(any("ENST00000999999" in h for h in headers))

    def test_sequence_body_preserved(self):
        filter_lncrna_fasta(self.fa, ["ENST00000761542"], self.out)
        with open(self.out) as f:
            body = f.read()
        self.assertIn("ACGTACGTAC\nACGTACGTAC", body)
        self.assertNotIn("TTTTTTTTTT", body)  # the other transcript's sequence

    def test_returns_output_path(self):
        ret = filter_lncrna_fasta(self.fa, ["ENST00000761542"], self.out)
        self.assertEqual(ret, self.out)


if __name__ == "__main__":
    unittest.main()
