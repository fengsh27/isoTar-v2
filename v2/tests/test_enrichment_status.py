"""Tests for app_v1.enrichment_status, the state file behind async enrichment.

Run with:
    python3 -m unittest v2.tests.test_enrichment_status
"""

import json
import os
import shutil
import sys
import tempfile
import unittest

# Make repo importable when running this file directly.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app_v1 import enrichment_status as es  # noqa: E402


def _touch(path, text=""):
    with open(path, "w") as f:
        f.write(text)


class EnrichmentStatusTests(unittest.TestCase):
    def setUp(self):
        self.job_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.job_dir)

    def test_never_run_reads_none(self):
        self.assertIsNone(es.read_status(self.job_dir))

    def test_legacy_run_reads_succeeded(self):
        # Enriched by the old synchronous API: gene list, no status file.
        _touch(os.path.join(self.job_dir, "gene_list.csv"), "gene_label\nTP53\n")
        state = es.read_status(self.job_dir)
        self.assertEqual(state["status"], "succeeded")
        self.assertTrue(state["legacy"])

    def test_write_then_update_round_trip(self):
        es.write_status(self.job_dir, {"status": "queued", "queued_at": 100,
                                       "params": {"organism": "Human"}})
        es.update_status(self.job_dir, status="running", started_at=200)
        state = es.read_status(self.job_dir)
        self.assertEqual(state["status"], "running")
        self.assertEqual(state["started_at"], 200)
        self.assertEqual(state["params"], {"organism": "Human"})  # kept by update
        self.assertIn("updated_at", state)
        self.assertFalse(os.path.exists(es.status_path(self.job_dir) + ".tmp"))

    def test_status_file_wins_over_legacy_gene_list(self):
        _touch(os.path.join(self.job_dir, "gene_list.csv"), "gene_label\nTP53\n")
        es.write_status(self.job_dir, {"status": "queued"})
        self.assertEqual(es.read_status(self.job_dir)["status"], "queued")

    def test_active_and_inactive(self):
        self.assertTrue(es.is_active({"status": "queued"}))
        self.assertTrue(es.is_active({"status": "running"}))
        self.assertFalse(es.is_active({"status": "succeeded"}))
        self.assertFalse(es.is_active({"status": "failed"}))
        self.assertFalse(es.is_active(None))

    def test_stale_running_reads_failed(self):
        es.write_status(self.job_dir, {"status": "running"})
        updated = es.read_status(self.job_dir)["updated_at"]
        fresh = es.read_status(self.job_dir, now=updated + 10)
        self.assertEqual(fresh["status"], "running")
        stale = es.read_status(self.job_dir, now=updated + es.STALE_AFTER + 1)
        self.assertEqual(stale["status"], "failed")
        self.assertIn("Run it again", stale["error"])
        self.assertFalse(es.is_active(stale))  # so a re-run is allowed
        # Reading never rewrites the file.
        with open(es.status_path(self.job_dir)) as f:
            self.assertEqual(json.load(f)["status"], "running")

    def test_finished_runs_never_go_stale(self):
        es.write_status(self.job_dir, {"status": "succeeded"})
        state = es.read_status(self.job_dir, now=10 ** 12)
        self.assertEqual(state["status"], "succeeded")

    def test_corrupt_file_reads_as_never_run(self):
        _touch(es.status_path(self.job_dir), "{not json")
        self.assertIsNone(es.read_status(self.job_dir))

    def test_clear_artifacts_removes_outputs_and_zip_only(self):
        keep = ["gene_list.csv", "job.json", "mirna.fa"]
        drop = ["KEGG_2021_Human.Human.enrichr.reports.txt",
                "GO_Biological_Process_2023.Human.enrichr.reports.txt",
                "enrichment_dotplot.png", "result.zip"]
        for name in keep + drop:
            _touch(os.path.join(self.job_dir, name))
        removed = sorted(os.path.basename(p) for p in es.clear_artifacts(self.job_dir))
        self.assertEqual(removed, sorted(drop))
        self.assertEqual(sorted(os.listdir(self.job_dir)), sorted(keep))


if __name__ == "__main__":
    unittest.main()
