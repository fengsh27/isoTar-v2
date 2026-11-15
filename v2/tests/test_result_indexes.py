"""Tests for app_v1/result_indexes.py: building a job's result databases in
the worker at job completion, and the warm-up CLI for existing jobs.

The ensure_* builders are stubbed, so these cover which indexes are planned
per workflow, error isolation, the rebuilt/current reporting and the CLI's
job discovery -- not the parsing itself (see test_result_db_sequences).

Needs the backend app deps (flask, celery, ...):
    python3 -m unittest v2.tests.test_result_indexes
"""

import io
import json
import os
import shutil
import sqlite3
import sys
import tempfile
import time
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app_v1 import result_indexes as ri  # noqa: E402


def _write_current_db(path):
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA user_version = {}".format(ri.SCHEMA_VERSION))
    conn.commit()
    conn.close()


def _write_old_db(path):
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE gene_tools (gene_id TEXT, tool TEXT)")
    conn.commit()
    conn.close()


class _FakeLogger(object):
    def __init__(self):
        self.records = []

    def info(self, msg, *args):
        self.records.append(("info", msg % args))

    def warning(self, msg, *args):
        self.records.append(("warning", msg % args))


class _StubbedBuilders(unittest.TestCase):
    """Replace the three ensure_* builders with recorders that create the file."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.calls = []
        self._orig = {n: getattr(ri, n) for n in
                      ("ensure_db", "ensure_lncrna_db", "ensure_network_db")}
        ri.ensure_db = self._builder("ensure_db", ri.DB_FILENAME)
        ri.ensure_lncrna_db = self._builder("ensure_lncrna_db", ri.LNCRNA_DB_FILENAME)
        ri.ensure_network_db = self._builder("ensure_network_db", ri.NETWORK_DB_FILENAME)

    def tearDown(self):
        for n, fn in self._orig.items():
            setattr(ri, n, fn)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _builder(self, name, filename):
        def build(directory):
            self.calls.append((name, directory))
            path = os.path.join(directory, filename)
            if not os.path.exists(path):
                _write_current_db(path)
            return path
        return build


class PlannedIndexesTest(_StubbedBuilders):

    def test_mir_target_builds_result_db_only(self):
        res = ri.build_result_indexes(self.tmp, "mir-target")
        self.assertEqual(self.calls, [("ensure_db", self.tmp)])
        self.assertEqual([r["name"] for r in res], ["result.db"])
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "result.db")))

    def test_mir_lncrna_builds_lncrna_db_only(self):
        ri.build_result_indexes(self.tmp, "mir-lncrna")
        self.assertEqual(self.calls, [("ensure_lncrna_db", self.tmp)])

    def test_unknown_workflow_falls_back_to_gene_pool(self):
        ri.build_result_indexes(self.tmp, None)
        self.assertEqual(self.calls, [("ensure_db", self.tmp)])

    def test_network_builds_both_pools_then_graph(self):
        gene = os.path.join(self.tmp, "gene")
        lnc = os.path.join(self.tmp, "lncrna")
        os.makedirs(gene)
        os.makedirs(lnc)
        res = ri.build_result_indexes(self.tmp, "mir-network")
        self.assertEqual(self.calls, [
            ("ensure_db", gene),
            ("ensure_lncrna_db", lnc),
            ("ensure_network_db", self.tmp),
        ])
        self.assertEqual([r["name"] for r in res],
                         ["gene/result.db", "lncrna/lncrna_result.db", "network.db"])

    def test_network_skips_a_pool_that_was_never_run(self):
        os.makedirs(os.path.join(self.tmp, "gene"))
        ri.build_result_indexes(self.tmp, "mir-network")
        self.assertEqual([c[0] for c in self.calls], ["ensure_db", "ensure_network_db"])


class ReportingTest(_StubbedBuilders):

    def test_reports_rebuilt_when_file_missing_or_old(self):
        _write_old_db(os.path.join(self.tmp, "result.db"))
        res = ri.build_result_indexes(self.tmp, "mir-target")
        self.assertTrue(res[0]["rebuilt"])
        self.assertIsNone(res[0]["error"])
        self.assertEqual(res[0]["path"], os.path.join(self.tmp, "result.db"))

    def test_reports_current_when_already_on_schema(self):
        _write_current_db(os.path.join(self.tmp, "result.db"))
        logger = _FakeLogger()
        res = ri.build_result_indexes(self.tmp, "mir-target", logger=logger, job_id="j1")
        self.assertFalse(res[0]["rebuilt"])
        # The builder is still called (it is the one that decides), cheaply.
        self.assertEqual(len(self.calls), 1)
        self.assertEqual(logger.records,
                         [("info", "result index current job_id=j1 index=result.db")])

    def test_network_db_is_current_whenever_it_exists(self):
        spec = ri.planned_indexes(self.tmp, "mir-network")[-1]
        self.assertEqual(spec["name"], "network.db")
        self.assertFalse(ri.is_current(spec))
        open(os.path.join(self.tmp, "network.db"), "w").close()
        self.assertTrue(ri.is_current(spec))

    def test_a_failing_index_is_reported_and_does_not_stop_the_others(self):
        os.makedirs(os.path.join(self.tmp, "gene"))
        os.makedirs(os.path.join(self.tmp, "lncrna"))

        def boom(directory):
            self.calls.append(("ensure_lncrna_db", directory))
            raise RuntimeError("no parameters json")
        ri.ensure_lncrna_db = boom

        logger = _FakeLogger()
        res = ri.build_result_indexes(self.tmp, "mir-network", logger=logger, job_id="j2")
        self.assertEqual([c[0] for c in self.calls],
                         ["ensure_db", "ensure_lncrna_db", "ensure_network_db"])
        self.assertEqual([r["error"] for r in res], [None, "no parameters json", None])
        self.assertFalse(res[1]["rebuilt"])
        self.assertTrue(res[0]["rebuilt"] and res[2]["rebuilt"])
        levels = [lvl for lvl, _ in logger.records]
        self.assertEqual(levels, ["info", "warning", "info"])
        self.assertIn("index=lncrna/lncrna_result.db", logger.records[1][1])

    def test_logs_build_time_for_a_rebuilt_index(self):
        logger = _FakeLogger()
        ri.build_result_indexes(self.tmp, "mir-target", logger=logger, job_id="j3")
        self.assertEqual(len(logger.records), 1)
        self.assertTrue(logger.records[0][1].startswith(
            "result index built job_id=j3 index=result.db seconds="))


class _JobsDir(_StubbedBuilders):
    """A fake ISOTAR_JOB_DIR with a few jobs in different states."""

    def setUp(self):
        super(_JobsDir, self).setUp()
        self.now = int(time.time())
        self._job("old-ok", "succeeded", "mir-target", self.now - 40 * 86400)
        self._job("new-ok", "succeeded", "mir-target", self.now - 3600)
        self._job("net-ok", "succeeded", "mir-network", self.now - 7200, pools=("gene", "lncrna"))
        self._job("failed", "failed", "mir-target", self.now - 100)
        self._job("running", "running", "mir-target", None)
        os.makedirs(os.path.join(self.tmp, "no-meta"))
        os.makedirs(os.path.join(self.tmp, "broken"))
        with open(os.path.join(self.tmp, "broken", "job.json"), "w") as fh:
            fh.write("{not json")
        # Succeeded but its output directory was deleted.
        self._job("gone", "succeeded", "mir-target", self.now - 10, output=False)

    def _job(self, job_id, status, workflow, finished_at, pools=(), output=True):
        job_dir = os.path.join(self.tmp, job_id)
        out = os.path.join(job_dir, "output")
        os.makedirs(job_dir)
        if output:
            os.makedirs(out)
            for p in pools:
                os.makedirs(os.path.join(out, p))
        meta = {"status": status, "workflow": workflow, "result_path": out}
        if finished_at is not None:
            meta["finished_at"] = finished_at
        with open(os.path.join(job_dir, "job.json"), "w") as fh:
            json.dump(meta, fh)
        return out


class DiscoverJobsTest(_JobsDir):

    def test_only_succeeded_jobs_with_output_newest_first(self):
        jobs = ri.discover_jobs(self.tmp)
        self.assertEqual([j["job_id"] for j in jobs], ["new-ok", "net-ok", "old-ok"])
        self.assertEqual(jobs[1]["workflow"], "mir-network")
        self.assertEqual(jobs[0]["output_dir"], os.path.join(self.tmp, "new-ok", "output"))

    def test_only_filter(self):
        jobs = ri.discover_jobs(self.tmp, only=["old-ok", "failed", "nope"])
        self.assertEqual([j["job_id"] for j in jobs], ["old-ok"])

    def test_since_filter(self):
        jobs = ri.discover_jobs(self.tmp, since_days=7)
        self.assertEqual([j["job_id"] for j in jobs], ["new-ok", "net-ok"])


class CliTest(_JobsDir):

    def _run(self, *argv):
        out = io.StringIO()
        rc = ri.main(["--jobs-dir", self.tmp] + list(argv), out=out)
        return rc, out.getvalue()

    def test_dry_run_lists_stale_indexes_and_builds_nothing(self):
        _write_current_db(os.path.join(self.tmp, "new-ok", "output", "result.db"))
        rc, text = self._run("--dry-run")
        self.assertEqual(rc, 0)
        self.assertEqual(self.calls, [])
        lines = text.strip().splitlines()
        self.assertEqual(lines[:4], [
            "stale    net-ok mir-network gene/result.db",
            "stale    net-ok mir-network lncrna/lncrna_result.db",
            "stale    net-ok mir-network network.db",
            "stale    old-ok mir-target result.db",
        ])
        self.assertEqual(lines[-1], "3 succeeded job(s), 4 stale index(es)")

    def test_build_reports_each_index_and_the_slowest(self):
        rc, text = self._run()
        self.assertEqual(rc, 0)
        self.assertEqual(len(self.calls), 5)
        lines = text.strip().splitlines()
        self.assertTrue(lines[0].startswith("built    new-ok mir-target result.db "))
        self.assertEqual(lines[-2], "3 succeeded job(s), 5 index(es) built, 0 failed")
        self.assertTrue(lines[-1].startswith("slowest build: "))

    def test_build_failure_sets_exit_code(self):
        def boom(directory):
            raise RuntimeError("disk full")
        ri.ensure_network_db = boom
        rc, text = self._run("--job", "net-ok")
        self.assertEqual(rc, 1)
        self.assertIn("FAILED   net-ok mir-network network.db", text)
        self.assertIn("disk full", text)
        self.assertIn("1 succeeded job(s), 2 index(es) built, 1 failed", text)

    def test_missing_jobs_dir(self):
        out = io.StringIO()
        rc = ri.main(["--jobs-dir", os.path.join(self.tmp, "nowhere")], out=out)
        self.assertEqual(rc, 2)
        self.assertIn("jobs dir not found", out.getvalue())


if __name__ == "__main__":
    unittest.main()
