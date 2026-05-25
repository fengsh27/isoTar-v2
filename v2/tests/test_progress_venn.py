"""Regression tests for two result-display fixes.

1. _write_progress (v2/mirna_predicting.py) must MERGE with the existing
   progress.json, because prediction runs in two phases (python3.6 for most
   tools, python2.7 for miRmap) that share the file. Before the fix the second
   phase clobbered the first, so a finished job's progress.json listed only
   miRmap and the UI showed the other tools as "pending" after a refresh.

2. _venn_stats (app_v1/result_db.py) must emit exclusive `combinations` and
   per-degree `degrees` so the 4+-tool UpSet plot and consensus histogram have
   data. Before the fix it returned only `sets` + `intersections`, so those
   views rendered empty / fell back to the flat table.

Runnable with:  python3 -m unittest v2.tests.test_progress_venn
"""
import json
import os
import shutil
import sqlite3
import sys
import tempfile
import unittest

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from v2.mirna_predicting import _write_progress
from app_v1.result_db import _venn_stats


def _done(start, end):
    return {"status": "done", "started_at": start, "finished_at": end}


class WriteProgressMergeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _read(self):
        with open(os.path.join(self.tmp, "progress.json")) as f:
            return json.load(f)

    def test_second_phase_does_not_clobber_first(self):
        # Phase 1 (python3.6): five tools finish.
        phase1 = {
            "miRanda": _done(1, 2),
            "Targetscan": _done(2, 3),
            "RNAhybrid": _done(3, 4),
            "PITA": _done(4, 5),
            "DMISO": _done(5, 6),
        }
        _write_progress(self.tmp, list(phase1.keys()), phase1)

        # Phase 2 (python2.7): miRmap only -- this used to overwrite the file.
        _write_progress(self.tmp, ["miRmap"], {"miRmap": _done(7, 8)})

        data = self._read()
        status = data["tools_status"]
        self.assertEqual(
            set(status.keys()),
            {"miRanda", "Targetscan", "RNAhybrid", "PITA", "DMISO", "miRmap"},
        )
        self.assertEqual(data["total_tools"], 6)
        self.assertEqual(data["completed_tools"], 6)
        # First-phase timing is preserved, not reset.
        self.assertEqual(status["miRanda"]["started_at"], 1)
        self.assertEqual(status["miRmap"]["finished_at"], 8)

    def test_running_tool_reported_as_current(self):
        _write_progress(self.tmp, ["miRanda"], {"miRanda": _done(1, 2)})
        _write_progress(
            self.tmp, ["miRmap"],
            {"miRmap": {"status": "running", "started_at": 3, "finished_at": None}},
        )
        data = self._read()
        self.assertEqual(data["current_tool"], "miRmap")
        self.assertEqual(data["completed_tools"], 1)
        self.assertEqual(data["total_tools"], 2)


class VennStatsTest(unittest.TestCase):
    def _cursor(self, rows):
        conn = sqlite3.connect(":memory:")
        c = conn.cursor()
        c.execute("CREATE TABLE gene_tools (gene_id TEXT, tool TEXT, "
                  "PRIMARY KEY (gene_id, tool))")
        c.executemany("INSERT INTO gene_tools VALUES (?, ?)", rows)
        c.execute("CREATE TABLE gene_info (gene_id TEXT PRIMARY KEY, "
                  "gene_label TEXT, gene_name TEXT)")
        conn.commit()
        return c

    def test_exclusive_combinations_and_degrees(self):
        # g1,g5 -> A,B,C ; g2 -> A,B ; g3 -> A ; g4 -> B,C
        rows = [
            ("g1", "A"), ("g1", "B"), ("g1", "C"),
            ("g2", "A"), ("g2", "B"),
            ("g3", "A"),
            ("g4", "B"), ("g4", "C"),
            ("g5", "A"), ("g5", "B"), ("g5", "C"),
        ]
        venn = _venn_stats(self._cursor(rows))

        self.assertEqual(venn["sets"], {"A": 4, "B": 4, "C": 3})

        combo = {tuple(c["tools"]): c["size"] for c in venn["combinations"]}
        self.assertEqual(combo[("A", "B", "C")], 2)
        self.assertEqual(combo[("A", "B")], 1)
        self.assertEqual(combo[("A",)], 1)
        self.assertEqual(combo[("B", "C")], 1)

        # Exclusive regions partition the gene set exactly once.
        self.assertEqual(sum(c["size"] for c in venn["combinations"]), 5)
        # Sorted by size, descending.
        sizes = [c["size"] for c in venn["combinations"]]
        self.assertEqual(sizes, sorted(sizes, reverse=True))

        self.assertEqual(venn["degrees"], {"1": 1, "2": 2, "3": 2})


if __name__ == "__main__":
    unittest.main()
