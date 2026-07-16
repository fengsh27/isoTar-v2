"""Tests for per-tool "time costed" accounting.

Every tool runs once per miRNA, and a mir-network job runs every tool over two
target pools (gene + lncRNA). The runner used to store a single started_at /
finished_at per tool, overwriting them on each miRNA, and the backend then
merged the pools with min(start)/max(finish). "Time costed" was read off that
bracket -- but a tool's bracket also spans every OTHER tool's work, so on a real
7h15m job (9b4b222b) the UI reported 3h58m for a miRanda that ran 66 seconds,
and the column summed to 21.4h. Duration now comes from `elapsed`, accumulated
per run by the runner and summed across pools by the backend.

Needs the backend app deps (flask, celery, ...), so run inside the backend
container (or any env with the app_v1 stack installed):
    python3 -m unittest v2.tests.test_progress_elapsed
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
from v2 import mirna_predicting as mp  # noqa: E402


class _FakeClock(object):
    """Stands in for the `time` module inside mirna_predicting."""

    def __init__(self, start=1000):
        self.now = start

    def time(self):
        return self.now

    def advance(self, seconds):
        self.now += seconds


class _ClockTestCase(unittest.TestCase):
    def setUp(self):
        self.clock = _FakeClock()
        self._real_time = mp.time
        mp.time = self.clock
        self.tmp = tempfile.mkdtemp(prefix="isotar_elapsed_")

    def tearDown(self):
        mp.time = self._real_time
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run_tool(self, statuses, tool, seconds):
        """One run of `tool` taking `seconds`."""
        mp._tool_started(statuses, tool)
        self.clock.advance(seconds)
        mp._tool_finished(statuses, tool)


class ToolElapsedAccumulationTest(_ClockTestCase):
    def test_elapsed_accumulates_across_miRNAs(self):
        statuses = mp._init_tool_statuses(["miRanda"])
        # Three miRNAs, 10s of work each, with a long unrelated gap between --
        # the gap is another tool's work and must not be attributed to miRanda.
        for _ in range(3):
            self._run_tool(statuses, "miRanda", 10)
            self.clock.advance(3600)

        info = statuses["miRanda"]
        self.assertEqual(info["elapsed"], 30)  # work, not span
        self.assertEqual(info["status"], "done")
        self.assertIsNone(info["running_since"])
        # The bracket is over two hours wide; that is exactly why it is not used.
        self.assertGreater(info["finished_at"] - info["started_at"], 7000)

    def test_started_at_keeps_the_first_start(self):
        statuses = mp._init_tool_statuses(["PITA"])
        self._run_tool(statuses, "PITA", 5)
        first_start = statuses["PITA"]["started_at"]
        self.clock.advance(500)
        self._run_tool(statuses, "PITA", 5)
        # Not overwritten by the later run -- "Started" means first start.
        self.assertEqual(statuses["PITA"]["started_at"], first_start)
        self.assertEqual(statuses["PITA"]["elapsed"], 10)

    def test_running_since_tracks_the_in_flight_run(self):
        statuses = mp._init_tool_statuses(["DMISO"])
        self._run_tool(statuses, "DMISO", 7)
        self.clock.advance(100)
        mp._tool_started(statuses, "DMISO")  # second run, still going
        info = statuses["DMISO"]
        self.assertEqual(info["status"], "running")
        self.assertEqual(info["elapsed"], 7)  # only the banked run
        self.assertEqual(info["running_since"], self.clock.now)

    def test_progress_json_carries_elapsed(self):
        statuses = mp._init_tool_statuses(["miRanda"])
        self._run_tool(statuses, "miRanda", 12)
        mp._write_progress(self.tmp, ["miRanda"], statuses)
        with open(os.path.join(self.tmp, "progress.json")) as f:
            data = json.load(f)
        self.assertEqual(data["tools_status"]["miRanda"]["elapsed"], 12)


class NetworkJobElapsedTest(_ClockTestCase):
    """Replays job 9b4b222b's real timings and asserts the UI would now show
    each tool's actual work instead of the overlapping bracket it did show."""

    # Seconds of real work per run, from the job's output-file mtimes.
    GENE = {
        "miRanda":   [11, 11, 11],
        "RNAhybrid": [3186, 366, 3066],
        "PITA":      [276, 264, 264],
        "DMISO":     [126, 126, 132],
    }
    LNCRNA = {
        "miRanda":   [11, 11, 11],
        "RNAhybrid": [4020, 414, 4002],
        "PITA":      [264, 264, 264],
        "DMISO":     [150, 156, 138],
    }
    GENE_MIRMAP = [1992, 2712, 834]
    LNCRNA_MIRMAP = [1134, 1122, 804]

    # What the frontend actually displayed for this job (the bug).
    OBSERVED_BRACKETS = {
        "miRanda": 14304,    # 3h 58m -- for 33s of real work
        "RNAhybrid": 18291,  # 5h 04m
        "PITA": 15490,       # 4h 18m
        "DMISO": 15366,      # 4h 16m
        "miRmap": 13594,     # 3h 46m
    }

    def _run_pool(self, pool_dir, per_tool, mirmap_runs):
        """Phase 1 (python3.6 tools) then phase 2 (miRmap), as the app does --
        two processes, two tool_statuses dicts, one shared progress.json."""
        os.makedirs(pool_dir)
        tools = ["miRanda", "RNAhybrid", "PITA", "DMISO"]
        statuses = mp._init_tool_statuses(tools)
        for i in range(3):  # three miRNAs
            for tool in tools:
                self._run_tool(statuses, tool, per_tool[tool][i])
                mp._write_progress(pool_dir, tools, statuses)

        mirmap_statuses = mp._init_tool_statuses(["miRmap"])
        for seconds in mirmap_runs:
            self._run_tool(mirmap_statuses, "miRmap", seconds)
            mp._write_progress(pool_dir, ["miRmap"], mirmap_statuses)

    def test_merged_elapsed_is_work_not_bracket(self):
        job_id = "net-elapsed"
        orig_base = app_module.BASE_DIR
        app_module.BASE_DIR = self.tmp
        try:
            output = os.path.join(app_module._job_path(job_id), "output")
            self._run_pool(os.path.join(output, "gene"), self.GENE, self.GENE_MIRMAP)
            self._run_pool(os.path.join(output, "lncrna"), self.LNCRNA, self.LNCRNA_MIRMAP)

            status = app_module._load_progress(job_id)["tools_status"]

            expected = {
                "miRanda": sum(self.GENE["miRanda"]) + sum(self.LNCRNA["miRanda"]),
                "RNAhybrid": sum(self.GENE["RNAhybrid"]) + sum(self.LNCRNA["RNAhybrid"]),
                "PITA": sum(self.GENE["PITA"]) + sum(self.LNCRNA["PITA"]),
                "DMISO": sum(self.GENE["DMISO"]) + sum(self.LNCRNA["DMISO"]),
                "miRmap": sum(self.GENE_MIRMAP) + sum(self.LNCRNA_MIRMAP),
            }
            for tool, want in expected.items():
                self.assertEqual(status[tool]["status"], "done")
                self.assertEqual(status[tool]["elapsed"], want)
                # And it is no longer the bracket the UI used to show.
                bracket = status[tool]["finished_at"] - status[tool]["started_at"]
                self.assertNotEqual(status[tool]["elapsed"], bracket)

            # miRanda: 66s of work, but its bracket spans ~4h of other tools'.
            self.assertEqual(status["miRanda"]["elapsed"], 66)
            self.assertGreater(
                status["miRanda"]["finished_at"] - status["miRanda"]["started_at"],
                10000)

            # The headline symptom: bracketed times summed to ~3x the job's
            # wall clock because they overlap. Real work cannot.
            total_work = sum(status[t]["elapsed"] for t in expected)
            wall = self.clock.now - 1000
            self.assertLessEqual(total_work, wall)
            self.assertGreater(sum(self.OBSERVED_BRACKETS.values()), wall * 2)
        finally:
            app_module.BASE_DIR = orig_base

    def test_matches_the_measured_totals(self):
        # Anchored to the per-tool totals measured from the job's own outputs.
        self.assertEqual(sum(self.GENE["RNAhybrid"]) + sum(self.LNCRNA["RNAhybrid"]), 15054)
        self.assertEqual(sum(self.GENE_MIRMAP) + sum(self.LNCRNA_MIRMAP), 8598)
        # 11s per run, 3 miRNAs, both pools -- against a displayed 3h 58m.
        self.assertEqual(sum(self.GENE["miRanda"]) + sum(self.LNCRNA["miRanda"]), 66)


class CombineElapsedTest(unittest.TestCase):
    """_combine_tool_status sums elapsed across pools, and stays honest about
    jobs that predate the field."""

    def _done(self, start, end, elapsed):
        return {"status": "done", "started_at": start, "finished_at": end,
                "elapsed": elapsed, "running_since": None}

    def test_sums_elapsed_across_pools(self):
        s = app_module._combine_tool_status(
            [self._done(1, 100, 30), self._done(200, 400, 45)])
        self.assertEqual(s["elapsed"], 75)          # work
        self.assertEqual(s["started_at"], 1)        # bracket start
        self.assertEqual(s["finished_at"], 400)     # bracket end
        self.assertNotEqual(s["elapsed"], s["finished_at"] - s["started_at"])

    def test_legacy_entries_have_no_elapsed(self):
        # Pre-`elapsed` progress.json: report None rather than inventing 0, so
        # the UI can fall back instead of claiming every tool took no time.
        legacy = {"status": "done", "started_at": 1, "finished_at": 100}
        s = app_module._combine_tool_status([legacy])
        self.assertIsNone(s["elapsed"])

    def test_running_since_is_earliest_in_flight_run(self):
        running = {"status": "running", "started_at": 5, "finished_at": None,
                   "elapsed": 10, "running_since": 50}
        s = app_module._combine_tool_status([self._done(1, 4, 3), running])
        self.assertEqual(s["status"], "running")
        self.assertEqual(s["elapsed"], 13)
        self.assertEqual(s["running_since"], 50)
        self.assertIsNone(s["finished_at"])


class FinalizeBanksElapsedTest(unittest.TestCase):
    """A tool whose runner died mid-run is flipped to 'failed'; its partial run
    must be banked and running_since cleared, or the UI live-counts a dead tool
    forever."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="isotar_finalize_")
        self._orig = app_module.BASE_DIR
        app_module.BASE_DIR = self.tmp

    def tearDown(self):
        app_module.BASE_DIR = self._orig
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_dead_run_is_banked_and_closed(self):
        job_id = "fin-elapsed"
        out = os.path.join(app_module._job_path(job_id), "output")
        os.makedirs(out)
        with open(os.path.join(out, "progress.json"), "w") as f:
            json.dump({"tools_status": {
                "Targetscan": {"status": "running", "started_at": 10,
                               "finished_at": None, "elapsed": 60,
                               "running_since": 100},
            }}, f)

        app_module._finalize_progress(job_id)

        with open(os.path.join(out, "progress.json")) as f:
            info = json.load(f)["tools_status"]["Targetscan"]
        self.assertEqual(info["status"], "failed")
        self.assertIsNone(info["running_since"])   # stops the live counter
        self.assertGreater(info["elapsed"], 60)    # partial run banked
        self.assertIsNotNone(info["finished_at"])


if __name__ == "__main__":
    unittest.main()
