"""run_job builds the result indexes in the worker before marking a job
succeeded, under an "indexing" step the status page can show.

Prediction itself is stubbed; this covers the step sequence and the hand-off
to app_v1.result_indexes, not the tools.

Needs the backend app deps (flask, celery, ...):
    python3 -m unittest v2.tests.test_run_job_indexing
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


class RunJobIndexingTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.job_id = "job-idx"
        os.makedirs(os.path.join(self.tmp, self.job_id))
        self._orig = {
            "BASE_DIR": app_module.BASE_DIR,
            "_process_single_mirna": app_module._process_single_mirna,
            "_process_mirna_list": app_module._process_mirna_list,
            "_run_prediction_pool": app_module._run_prediction_pool,
            "_finalize_progress": app_module._finalize_progress,
            "build_result_indexes": app_module.build_result_indexes,
        }
        app_module.BASE_DIR = self.tmp
        app_module._process_single_mirna = lambda *a, **k: None
        app_module._process_mirna_list = lambda *a, **k: None
        app_module._finalize_progress = lambda job_id: None
        self.pool_calls = []
        self.index_calls = []

        def fake_pool(job_id, fasta_path, tools, genome, target_type, out_dir,
                      cores, target_file, tool_errors, label=None):
            self.pool_calls.append((target_type, out_dir, label))
            os.makedirs(out_dir, exist_ok=True)
            return 1
        app_module._run_prediction_pool = fake_pool

        def fake_build(output_dir, workflow, logger=None, job_id=None):
            # Record what the status endpoint would have shown at this moment.
            with open(os.path.join(self.tmp, job_id, "job.json")) as fh:
                meta = json.load(fh)
            self.index_calls.append((output_dir, workflow, meta["status"], meta["step"]))
            return []
        app_module.build_result_indexes = fake_build

    def tearDown(self):
        for name, value in self._orig.items():
            setattr(app_module, name, value)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write_meta(self, **extra):
        meta = {"status": "queued", "tools": ["miRanda"], "genome": "hg38",
                "mirna_id": "hsa-miR-21-5p", "cores": 1}
        meta.update(extra)
        with open(os.path.join(self.tmp, self.job_id, "job.json"), "w") as fh:
            json.dump(meta, fh)

    def _read_meta(self):
        with open(os.path.join(self.tmp, self.job_id, "job.json")) as fh:
            return json.load(fh)

    def test_indexes_built_while_running_then_succeeded(self):
        self._write_meta(workflow="mir-target")
        app_module.run_job(self.job_id)
        output_dir = os.path.join(self.tmp, self.job_id, "output")
        self.assertEqual(self.index_calls,
                         [(output_dir, "mir-target", "running", "indexing")])
        meta = self._read_meta()
        self.assertEqual(meta["status"], "succeeded")
        self.assertEqual(meta["result_path"], output_dir)
        self.assertIn("finished_at", meta)

    def test_network_job_passes_its_workflow(self):
        self._write_meta(workflow="mir-network", mirna_ids=["hsa-miR-21-5p"])
        app_module.run_job(self.job_id)
        self.assertEqual([c[1] for c in self.index_calls], ["mir-network"])
        self.assertEqual([c[2] for c in self.pool_calls], ["gene", "lncrna"])
        self.assertEqual(self._read_meta()["status"], "succeeded")

    def test_no_indexing_when_every_tool_group_failed(self):
        def failing_pool(job_id, fasta_path, tools, genome, target_type, out_dir,
                         cores, target_file, tool_errors, label=None):
            tool_errors[label or "pool"] = "boom"
            return 1
        app_module._run_prediction_pool = failing_pool
        self._write_meta(workflow="mir-target")
        app_module.run_job(self.job_id)
        self.assertEqual(self.index_calls, [])
        self.assertEqual(self._read_meta()["status"], "failed")


if __name__ == "__main__":
    unittest.main()
