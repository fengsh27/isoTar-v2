"""Tests for the `sequence` field of POST /api/v1/jobs/<id>/enrichment.

When the task collects the job's own genes it reads result.db, which keeps
WT and variant hits apart; `sequence` names the record to enrich and defaults
to the WT. The task itself is not run here (it shells out to the enrichment
script); the route's validation and the recorded params are.

Needs the backend app deps (flask, celery), so run inside the backend
container (or any env with the app_v1 stack installed):
    python3 -m unittest v2.tests.test_enrichment_route_sequence
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

try:
    import flask  # noqa: F401
    import celery  # noqa: F401
    _HAS_WEB_DEPS = True
except ImportError:
    _HAS_WEB_DEPS = False


class _FakeTask(object):
    def __init__(self):
        self.queued = []

    def apply_async(self, args=None, queue=None):
        self.queued.append((args, queue))
        return type("Result", (), {"id": "task-1"})()


@unittest.skipUnless(_HAS_WEB_DEPS, "flask/celery not installed outside the image")
class EnrichmentSequenceTests(unittest.TestCase):
    def setUp(self):
        self.base = tempfile.mkdtemp(prefix="isotar-enrich-")
        from app_v1 import app as appmod
        self.appmod = appmod
        self._saved_base = appmod.BASE_DIR
        self._saved_task = appmod.run_enrichment_task
        appmod.BASE_DIR = self.base
        self.task = appmod.run_enrichment_task = _FakeTask()
        appmod.app.config["TESTING"] = True
        self.client = appmod.app.test_client()

        job_dir = os.path.join(self.base, "j1")
        out_dir = os.path.join(job_dir, "output")
        os.makedirs(out_dir)
        with open(os.path.join(job_dir, "job.json"), "w") as f:
            json.dump({"job_id": "j1", "status": "succeeded", "workflow": "mir-target",
                       "result_path": out_dir}, f)

    def tearDown(self):
        self.appmod.BASE_DIR = self._saved_base
        self.appmod.run_enrichment_task = self._saved_task
        shutil.rmtree(self.base, ignore_errors=True)

    def _post(self, body):
        return self.client.post("/api/v1/jobs/j1/enrichment", data=json.dumps(body),
                                content_type="application/json")

    def test_default_records_no_sequence(self):
        resp = self._post({})
        self.assertEqual(resp.status_code, 202, resp.data)
        params = json.loads(resp.data.decode("utf-8"))["params"]
        self.assertEqual(params["gene_source"], "job")
        self.assertIsNone(params["sequence"])
        self.assertEqual(self.task.queued[0][0][1]["sequence"], None)

    def test_sequence_is_recorded_for_the_task(self):
        resp = self._post({"sequence": " hsa-miR-21-5p,1|0,shifted "})
        self.assertEqual(resp.status_code, 202, resp.data)
        params = json.loads(resp.data.decode("utf-8"))["params"]
        self.assertEqual(params["sequence"], "hsa-miR-21-5p,1|0,shifted")
        self.assertEqual(self.task.queued[0][0][1]["sequence"], "hsa-miR-21-5p,1|0,shifted")

    def test_sequence_must_be_a_non_empty_string(self):
        self.assertEqual(self._post({"sequence": 5}).status_code, 400)
        self.assertEqual(self._post({"sequence": "  "}).status_code, 400)
        self.assertEqual(self.task.queued, [])

    def test_sequence_is_rejected_with_an_explicit_gene_list(self):
        resp = self._post({"genes": ["TP53"], "sequence": "hsa-miR-21-5p"})
        self.assertEqual(resp.status_code, 400)
        self.assertEqual(self.task.queued, [])


if __name__ == "__main__":
    unittest.main()
