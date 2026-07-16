"""Tests for /api/v1/jobs/<id>/result across workflows.

A mir-network job predicts against two target pools (gene + lncRNA), each in its
own output subdirectory with its own mirna_prediction_parameters.json. The route
used to fall through to the single-pool path and read that file from the
top-level output dir, where a network job never writes one -- so /result raised
FileNotFoundError and returned 500 for every network job. It now takes a `pool`
argument and reads the matching subdirectory.

Needs the backend app deps (flask, celery, ...), so run inside the backend
container (or any env with the app_v1 stack installed):
    python3 -m unittest v2.tests.test_result_pools
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


_FAKE_PAGE = {
    "total_genes": 2,
    "total": 2,
    "genes": [],
    "venn": {"sets": {}, "intersections": {}},
}


class _ResultRouteTestCase(unittest.TestCase):
    """Stubs the db builders so these tests cover routing, not parsing."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="isotar_result_pools_")
        self._orig_base = app_module.BASE_DIR
        app_module.BASE_DIR = self.tmp

        self.calls = []
        self._orig = {
            "ensure_db": app_module.ensure_db,
            "ensure_lncrna_db": app_module.ensure_lncrna_db,
            "query_genes": app_module.query_genes,
        }

        def fake_ensure_db(output_dir):
            self.calls.append(("ensure_db", output_dir))
            return os.path.join(output_dir, "result.db")

        def fake_ensure_lncrna_db(output_dir):
            self.calls.append(("ensure_lncrna_db", output_dir))
            return os.path.join(output_dir, "lncrna_result.db")

        def fake_query_genes(db_path, **kwargs):
            self.calls.append(("query_genes", db_path))
            return dict(_FAKE_PAGE)

        app_module.ensure_db = fake_ensure_db
        app_module.ensure_lncrna_db = fake_ensure_lncrna_db
        app_module.query_genes = fake_query_genes

        app_module.app.config["TESTING"] = True
        self.client = app_module.app.test_client()

    def tearDown(self):
        app_module.BASE_DIR = self._orig_base
        for name, fn in self._orig.items():
            setattr(app_module, name, fn)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _make_job(self, job_id, workflow, pools=()):
        job_dir = app_module._job_path(job_id)
        output_dir = os.path.join(job_dir, "output")
        os.makedirs(output_dir)
        for pool in pools:
            os.makedirs(os.path.join(output_dir, pool))
        with open(app_module._job_meta_path(job_id), "w") as f:
            json.dump({
                "job_id": job_id,
                "status": "succeeded",
                "workflow": workflow,
                "result_path": output_dir,
            }, f)
        return output_dir

    def _get(self, job_id, query=""):
        return self.client.get("/api/v1/jobs/{}/result{}".format(job_id, query))


class NetworkResultPoolTest(_ResultRouteTestCase):
    def test_missing_pool_is_a_clear_400_not_a_500(self):
        # The reported bug: this returned 500 "JSON file not found:
        # .../output/mirna_prediction_parameters.json".
        self._make_job("net", "mir-network", pools=("gene", "lncrna"))
        resp = self._get("net")
        self.assertEqual(resp.status_code, 400)
        body = json.loads(resp.data.decode("utf-8"))
        self.assertIn("pool", body["error"])
        self.assertEqual(body["pools"], ["gene", "lncrna"])
        self.assertIn("/network", body["hint"])  # points at the graph endpoint
        # Nothing was parsed.
        self.assertEqual(self.calls, [])

    def test_invalid_pool_rejected(self):
        self._make_job("net", "mir-network", pools=("gene", "lncrna"))
        resp = self._get("net", "?pool=mirna")
        self.assertEqual(resp.status_code, 400)
        self.assertEqual(self.calls, [])

    def test_gene_pool_reads_the_gene_subdirectory(self):
        output_dir = self._make_job("net", "mir-network", pools=("gene", "lncrna"))
        resp = self._get("net", "?pool=gene")
        self.assertEqual(resp.status_code, 200)
        body = json.loads(resp.data.decode("utf-8"))
        self.assertEqual(body["pool"], "gene")
        self.assertEqual(body["total_genes"], 2)
        # The gene pool is shaped like mir-target -> the RefSeq/symbol builder,
        # pointed at output/gene rather than output/.
        self.assertEqual(self.calls[0],
                         ("ensure_db", os.path.join(output_dir, "gene")))

    def test_lncrna_pool_uses_the_transcript_builder(self):
        output_dir = self._make_job("net", "mir-network", pools=("gene", "lncrna"))
        resp = self._get("net", "?pool=lncrna")
        self.assertEqual(resp.status_code, 200)
        self.assertEqual(json.loads(resp.data.decode("utf-8"))["pool"], "lncrna")
        self.assertEqual(self.calls[0],
                         ("ensure_lncrna_db", os.path.join(output_dir, "lncrna")))

    def test_absent_pool_directory_is_404(self):
        # Gene-only network run: asking for lncrna is missing data, not a crash.
        self._make_job("net", "mir-network", pools=("gene",))
        resp = self._get("net", "?pool=lncrna")
        self.assertEqual(resp.status_code, 404)
        self.assertEqual(self.calls, [])

    def test_pool_survives_the_other_query_params(self):
        self._make_job("net", "mir-network", pools=("gene", "lncrna"))
        resp = self._get("net", "?pool=gene&sortBy=gene_label&number=5&offset=10")
        self.assertEqual(resp.status_code, 200)
        body = json.loads(resp.data.decode("utf-8"))
        self.assertEqual(body["pool"], "gene")
        self.assertEqual(body["sort_by"], "gene_label")
        self.assertEqual(body["number"], 5)
        self.assertEqual(body["offset"], 10)


class SinglePoolWorkflowsUnchangedTest(_ResultRouteTestCase):
    """The single-pool workflows must keep reading the top-level output dir."""

    def test_mir_target_uses_gene_builder_at_top_level(self):
        output_dir = self._make_job("tgt", "mir-target")
        resp = self._get("tgt")
        self.assertEqual(resp.status_code, 200)
        self.assertEqual(self.calls[0], ("ensure_db", output_dir))
        self.assertNotIn("pool", json.loads(resp.data.decode("utf-8")))

    def test_mir_lncrna_uses_transcript_builder_at_top_level(self):
        output_dir = self._make_job("lnc", "mir-lncrna")
        resp = self._get("lnc")
        self.assertEqual(resp.status_code, 200)
        self.assertEqual(self.calls[0], ("ensure_lncrna_db", output_dir))

    def test_pool_is_rejected_for_non_network_jobs(self):
        self._make_job("tgt", "mir-target")
        resp = self._get("tgt", "?pool=gene")
        self.assertEqual(resp.status_code, 400)
        self.assertEqual(self.calls, [])


class ResultPreconditionsTest(_ResultRouteTestCase):
    def test_unknown_job_is_404(self):
        self.assertEqual(self._get("nope").status_code, 404)

    def test_unfinished_job_is_409_before_pool_validation(self):
        job_dir = app_module._job_path("run")
        os.makedirs(os.path.join(job_dir, "output"))
        with open(app_module._job_meta_path("run"), "w") as f:
            json.dump({"job_id": "run", "status": "running",
                       "workflow": "mir-network",
                       "result_path": os.path.join(job_dir, "output")}, f)
        # Still 409 even though `pool` is absent -- job state outranks it.
        self.assertEqual(self._get("run").status_code, 409)


if __name__ == "__main__":
    unittest.main()
