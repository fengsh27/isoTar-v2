"""Tests for v2.mirna_predicting.parse_fasta's record typing.

The runner tags each FASTA record with the variant type named in its header.
"modified_shifted" contains both "modified" and "shifted", so it has to be
tested before either, or a combined variant is typed as plain "modified".

Run with:
    python3 -m unittest v2.tests.test_parse_fasta_types
"""

import os
import shutil
import sys
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from v2.mirna_predicting import parse_fasta  # noqa: E402


class ParseFastaTypeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="isotar_fasta_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _parse(self, records):
        path = os.path.join(self.tmp, "mirna.fa")
        with open(path, "w") as fh:
            for header, seq in records:
                fh.write(">{}\n{}\n".format(header, seq))
        return parse_fasta(path)

    def test_every_variant_type(self):
        seq = "UAGCUUAUCAGACUGAUGUUGA"
        recs = self._parse([
            ("hsa-miR-21-5p,WT", seq),
            ("hsa-miR-21-5p,8:A|U,modified", seq),
            ("hsa-miR-21-5p,-7|1,shifted", seq),
            ("hsa-miR-21-5p,8:A|U,-7|1,modified_shifted", seq),
            ("hsa-miR-21-5p", seq),
        ])
        self.assertEqual([r["type"] for r in recs],
                         ["WT", "modified", "shifted", "modified_shifted", "unknown"])
        self.assertEqual(recs[0]["header"], "hsa-miR-21-5p,WT")
        self.assertEqual(recs[0]["length"], len(seq))

    def test_length_bounds(self):
        recs = self._parse([
            ("short,WT", "ACGU" * 4),        # 16 nt: dropped
            ("ok,WT", "ACGU" * 5),           # 20 nt
            ("long,WT", "ACGU" * 8),         # 32 nt: dropped
        ])
        self.assertEqual([r["header"] for r in recs], ["ok,WT"])


if __name__ == "__main__":
    unittest.main()
