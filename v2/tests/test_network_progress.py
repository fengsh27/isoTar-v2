"""Tests for mir-network job progress aggregation in app_v1.app.

A network job runs every tool against two target pools (gene + lncRNA), each in
its own output subdirectory, so each pool writes its OWN progress.json and there
is no top-level output/progress.json. _load_progress must merge the per-pool
files into one per-tool view keyed by bare tool name; before the fix it read
only the (never-created) top-level file, so the UI showed every tool "pending".

Needs the backend app deps (flask, celery, ...), so run inside the backend
container (or any env with the app_v1 stack installed):
    python3 -m unittest v2.tests.test_network_progress
"""

import json
import os
import shutil
import sys
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app_v1 import app as app_module  # noqa: E402


def _done(start, end):
    return {"status": "done", "started_at": start, "finished_at": end}


def _running(start):
    return {"status": "running", "started_at": start, "finished_at": None}


def _pending():
    return {"status": "pending", "started_at": None, "finished_at": None}


class _BaseDirTestCase(unittest.TestCase):
    """Redirect job state into a temp dir so tests never touch /opt/out/jobs."""

    def setUp(self):
        self._tmp = tempfile.mkdtemp(prefix="isotar_test_jobs_")
        self._orig_base_dir = app_module.BASE_DIR
        app_module.BASE_DIR = self._tmp

    def tearDown(self):
        app_module.BASE_DIR = self._orig_base_dir
        shutil.rmtree(self._tmp, ignore_errors=True)

    def _write_pool(self, job_id, pool, tools_status):
        pool_dir = os.path.join(app_module._job_path(job_id), "output", pool)
        os.makedirs(pool_dir)
        with open(os.path.join(pool_dir, "progress.json"), "w") as f:
            json.dump({"tools_status": tools_status, "updated_at": 1}, f)

    def _write_top(self, job_id, payload):
        out_dir = os.path.join(app_module._job_path(job_id), "output")
        os.makedirs(out_dir)
        with open(os.path.join(out_dir, "progress.json"), "w") as f:
            json.dump(payload, f)


class LoadProgressNetworkMergeTest(_BaseDirTestCase):
    def test_merges_both_pools_keyed_by_bare_tool_name(self):
        job_id = "net-1"
        # gene pool ran all three tools; lncRNA pool ran the two compatible ones
        # (TargetScan is gene-only) and is still working through miRanda.
        self._write_pool(job_id, "gene", {
            "miRanda": _done(1, 4),
            "miRmap": _done(4, 6),
            "Targetscan": _done(6, 8),
        })
        self._write_pool(job_id, "lncrna", {
            "miRanda": _running(9),
            "miRmap": _pending(),
        })

        progress = app_module._load_progress(job_id)
        self.assertIsNotNone(progress)
        status = progress["tools_status"]
        # Keys are bare tool names (what the frontend matches job['tools'] on).
        self.assertEqual(set(status.keys()), {"miRanda", "miRmap", "Targetscan"})
        # miRanda is running in the lncRNA pool -> overall running.
        self.assertEqual(status["miRanda"]["status"], "running")
        self.assertEqual(status["miRanda"]["started_at"], 1)  # earliest start
        self.assertIsNone(status["miRanda"]["finished_at"])
        # miRmap done in gene but still pending in lncRNA -> not done yet.
        self.assertEqual(status["miRmap"]["status"], "pending")
        # TargetScan only ran in the gene pool and finished -> done.
        self.assertEqual(status["Targetscan"]["status"], "done")
        self.assertEqual(status["Targetscan"]["finished_at"], 8)

        self.assertEqual(progress["total_tools"], 3)
        self.assertEqual(progress["completed_tools"], 1)  # only TargetScan
        self.assertEqual(progress["current_tool"], "miRanda")

    def test_all_pools_done_marks_tool_done_with_spanning_timing(self):
        job_id = "net-2"
        self._write_pool(job_id, "gene", {"miRanda": _done(1, 5)})
        self._write_pool(job_id, "lncrna", {"miRanda": _done(6, 12)})

        status = app_module._load_progress(job_id)["tools_status"]
        self.assertEqual(status["miRanda"]["status"], "done")
        self.assertEqual(status["miRanda"]["started_at"], 1)   # earliest start
        self.assertEqual(status["miRanda"]["finished_at"], 12)  # latest finish

    def test_top_level_file_takes_precedence(self):
        # A non-network job has a real top-level progress.json; merge must not
        # kick in and shadow it.
        job_id = "single-1"
        self._write_top(job_id, {
            "total_tools": 1,
            "completed_tools": 1,
            "current_tool": None,
            "tools_status": {"miRanda": _done(1, 2)},
            "updated_at": 9,
        })
        progress = app_module._load_progress(job_id)
        self.assertEqual(progress["updated_at"], 9)
        self.assertEqual(progress["tools_status"]["miRanda"]["status"], "done")

    def test_no_progress_files_returns_none(self):
        job_id = "empty-1"
        os.makedirs(os.path.join(app_module._job_path(job_id), "output"))
        self.assertIsNone(app_module._load_progress(job_id))


class CombineToolStatusTest(unittest.TestCase):
    """_combine_tool_status folds a tool's per-pool states into one."""

    def test_running_wins(self):
        s = app_module._combine_tool_status([_done(1, 2), _running(3)])
        self.assertEqual(s["status"], "running")

    def test_all_done_is_done(self):
        s = app_module._combine_tool_status([_done(1, 5), _done(6, 9)])
        self.assertEqual(s["status"], "done")
        self.assertEqual(s["finished_at"], 9)

    def test_failed_beats_pending_when_not_running(self):
        failed = {"status": "failed", "started_at": 3, "finished_at": 4}
        s = app_module._combine_tool_status([failed, _pending()])
        self.assertEqual(s["status"], "failed")

    def test_pending_when_no_terminal_signal(self):
        s = app_module._combine_tool_status([_pending(), _pending()])
        self.assertEqual(s["status"], "pending")


class FinalizeProgressTest(_BaseDirTestCase):
    """Once every prediction group returns, no tool may stay 'running'/'pending':
    such a tool had its runner die mid-run, so it must show as 'failed'. This is
    the fix for a finished job displaying a tool stuck on 'Running'."""

    def _read_top(self, job_id):
        path = os.path.join(app_module._job_path(job_id), "output", "progress.json")
        with open(path) as f:
            return json.load(f)

    def test_top_level_running_is_flipped_to_failed(self):
        job_id = "fin-1"
        # Mirrors the observed bug: TargetScan crashed after being marked running
        # while the other tools finished; the job still succeeded on their output.
        self._write_top(job_id, {
            "total_tools": 6,
            "completed_tools": 5,
            "current_tool": None,
            "tools_status": {
                "miRanda": _done(1, 2),
                "Targetscan": _running(3),
                "RNAhybrid": _done(4, 5),
                "PITA": _done(5, 6),
                "DMISO": _done(6, 7),
                "miRmap": _done(8, 9),
            },
            "updated_at": 9,
        })
        app_module._finalize_progress(job_id)
        data = self._read_top(job_id)
        ts = data["tools_status"]
        self.assertEqual(ts["Targetscan"]["status"], "failed")
        self.assertIsNotNone(ts["Targetscan"]["finished_at"])
        # Every other tool is untouched.
        for tool in ("miRanda", "RNAhybrid", "PITA", "DMISO", "miRmap"):
            self.assertEqual(ts[tool]["status"], "done")
        self.assertIsNone(data["current_tool"])
        self.assertEqual(data["completed_tools"], 5)

    def test_pending_is_also_failed(self):
        job_id = "fin-2"
        self._write_top(job_id, {
            "tools_status": {"miRanda": _done(1, 2), "PITA": _pending()},
        })
        app_module._finalize_progress(job_id)
        ts = self._read_top(job_id)["tools_status"]
        self.assertEqual(ts["PITA"]["status"], "failed")

    def test_all_done_left_unchanged(self):
        job_id = "fin-3"
        payload = {
            "total_tools": 1,
            "completed_tools": 1,
            "current_tool": None,
            "tools_status": {"miRanda": _done(1, 2)},
            "updated_at": 7,
        }
        self._write_top(job_id, payload)
        app_module._finalize_progress(job_id)
        # Untouched: still done, updated_at not bumped.
        self.assertEqual(self._read_top(job_id), payload)

    def test_reconciles_each_network_pool_file(self):
        job_id = "fin-net"
        self._write_pool(job_id, "gene", {
            "miRanda": _done(1, 4),
            "Targetscan": _running(6),
        })
        self._write_pool(job_id, "lncrna", {
            "miRanda": _running(9),
        })
        app_module._finalize_progress(job_id)
        # Both pool files reconciled; the merged view then reports failures.
        status = app_module._load_progress(job_id)["tools_status"]
        self.assertEqual(status["Targetscan"]["status"], "failed")
        self.assertEqual(status["miRanda"]["status"], "failed")


if __name__ == "__main__":
    unittest.main()
