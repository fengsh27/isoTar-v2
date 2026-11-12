"""Tests for app_v1.result_db.gene_labels, the gene list enrichment runs on
when the client does not send one.

Run with:
    python3 -m unittest v2.tests.test_result_db_gene_labels
"""

import os
import shutil
import sqlite3
import sys
import tempfile
import unittest

# Make repo importable when running this file directly.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app_v1.result_db import gene_labels  # noqa: E402


class GeneLabelsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.db = os.path.join(self.tmp, "result.db")
        conn = sqlite3.connect(self.db)
        conn.executescript(
            """
            CREATE TABLE gene_tools (gene_id TEXT, tool TEXT, PRIMARY KEY (gene_id, tool));
            CREATE TABLE gene_info  (gene_id TEXT PRIMARY KEY, gene_label TEXT, gene_name TEXT);
            """)
        self.conn = conn

    def tearDown(self):
        self.conn.close()
        shutil.rmtree(self.tmp)

    def _add(self, gene_id, label, tools=("miRanda",)):
        for tool in tools:
            self.conn.execute("INSERT INTO gene_tools VALUES (?, ?)", (gene_id, tool))
        if label is not False:
            self.conn.execute("INSERT INTO gene_info VALUES (?, ?, ?)", (gene_id, label, ""))
        self.conn.commit()

    def test_transcripts_collapse_to_one_symbol(self):
        # Two RefSeq transcripts of TP53, each hit by several tools.
        self._add("NM_000546", "TP53", tools=("miRanda", "PITA"))
        self._add("NM_001126112", "TP53", tools=("RNAhybrid",))
        self._add("NM_002467", "MYC")
        self.assertEqual(gene_labels(self.db), ["MYC", "TP53"])

    def test_targets_without_a_symbol_are_skipped(self):
        self._add("NM_000546", "TP53")
        self._add("NM_999998", None)       # gene_info row with no label
        self._add("NM_999997", "  ")       # blank label
        self._add("NM_999996", False)      # no gene_info row at all
        self.assertEqual(gene_labels(self.db), ["TP53"])

    def test_labels_are_trimmed(self):
        self._add("NM_000546", " TP53 ")
        self.assertEqual(gene_labels(self.db), ["TP53"])

    def test_empty_result(self):
        self.assertEqual(gene_labels(self.db), [])


if __name__ == "__main__":
    unittest.main()
