"""Tests for the result download archive (GET /api/v1/jobs/<id>/result/download).

Covers what ends up inside result.zip, and in particular the cache
invalidation: the archive is built once and reused, but enrichment runs AFTER
the job succeeds, so a zip built before enrichment must not be served
afterwards.

Importing app_v1.app pulls in flask + celery, which are only present in the
container image -- the whole module skips when they are missing.

Run with:
    python3 -m unittest v2.tests.test_result_download_zip
"""

import json
import os
import shutil
import sys
import tempfile
import unittest
import zipfile

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


@unittest.skipUnless(_HAS_WEB_DEPS, "flask/celery not installed outside the image")
class ResultDownloadZipTests(unittest.TestCase):

    def setUp(self):
        self.base = tempfile.mkdtemp(prefix="isotar-jobs-")
        os.environ["ISOTAR_JOB_DIR"] = self.base
        # BASE_DIR is read at import time, so re-read it here rather than
        # relying on import order across tests.
        from app_v1 import app as appmod
        self.appmod = appmod
        self._saved_base = appmod.BASE_DIR
        appmod.BASE_DIR = self.base
        self.client = appmod.app.test_client()

    def tearDown(self):
        self.appmod.BASE_DIR = self._saved_base
        shutil.rmtree(self.base, ignore_errors=True)

    # ---- helpers ----

    def _make_job(self, job_id="j1"):
        job_dir = os.path.join(self.base, job_id)
        out_dir = os.path.join(job_dir, "output")
        os.makedirs(out_dir)
        with open(os.path.join(out_dir, "miranda.txt"), "w") as f:
            f.write("raw miranda output\n")
        with open(os.path.join(job_dir, "targets.txt"), "w") as f:
            f.write("NM_000546\n")
        with open(os.path.join(job_dir, "mirna.fa"), "w") as f:
            f.write(">hsa-miR-21-5p\nUAGCUUAUCAGACUGAUGUUGA\n")
        with open(os.path.join(job_dir, "job.json"), "w") as f:
            json.dump({"job_id": job_id, "status": "succeeded",
                       "result_path": out_dir}, f)
        return job_dir

    def _add_enrichment(self, job_dir, genes=("TP53",)):
        with open(os.path.join(job_dir, "gene_list.csv"), "w") as f:
            f.write("gene_label\n" + "".join(g + "\n" for g in genes))
        with open(os.path.join(job_dir,
                               "GO_Biological_Process_2023.Human.enrichr.reports.txt"), "w") as f:
            f.write("Gene_set\tTerm\tP-value\n")
        with open(os.path.join(job_dir, "enrichment_dotplot.png"), "wb") as f:
            f.write(b"\x89PNG\r\n\x1a\n")

    def _touch_newer(self, job_dir):
        """Push the enrichment artifacts past the zip's mtime without sleeping."""
        zip_mtime = os.path.getmtime(os.path.join(job_dir, "result.zip"))
        for name in os.listdir(job_dir):
            if name == "result.zip":
                continue
            os.utime(os.path.join(job_dir, name), (zip_mtime + 10, zip_mtime + 10))

    def _download(self, job_id="j1"):
        resp = self.client.get("/api/v1/jobs/{}/result/download".format(job_id))
        self.assertEqual(resp.status_code, 200, resp.data[:400])
        path = os.path.join(self.base, job_id + ".dl.zip")
        with open(path, "wb") as f:
            f.write(resp.data)
        with zipfile.ZipFile(path) as zf:
            return sorted(zf.namelist())

    # ---- tests ----

    def test_includes_results_inputs_and_enrichment(self):
        job_dir = self._make_job()
        self._add_enrichment(job_dir)
        names = self._download()
        self.assertIn("miranda.txt", names)  # result_dir/* at the top level
        for run_input in ("targets.txt", "job.json", "mirna.fa"):
            self.assertIn(run_input, names)
        self.assertEqual(
            [n for n in names if n.startswith("enrichment/")],
            ["enrichment/GO_Biological_Process_2023.Human.enrichr.reports.txt",
             "enrichment/enrichment_dotplot.png",
             "enrichment/gene_list.csv"])

    def test_no_enrichment_directory_when_not_run(self):
        self._make_job()
        names = self._download()
        self.assertFalse([n for n in names if n.startswith("enrichment/")])

    def test_cached_zip_rebuilt_when_enrichment_runs_later(self):
        job_dir = self._make_job()
        self.assertFalse([n for n in self._download() if n.startswith("enrichment/")])
        self._add_enrichment(job_dir)
        self._touch_newer(job_dir)
        self.assertTrue([n for n in self._download() if n.startswith("enrichment/")])

    def test_rerunning_enrichment_invalidates_the_cache(self):
        job_dir = self._make_job()
        self._add_enrichment(job_dir)
        self._download()
        zip_path = os.path.join(job_dir, "result.zip")
        before = os.path.getmtime(zip_path)
        self._add_enrichment(job_dir, genes=("TP53", "EGFR"))
        self._touch_newer(job_dir)
        self._download()
        self.assertGreater(os.path.getmtime(zip_path), before)

    def test_cache_reused_when_there_is_no_enrichment(self):
        job_dir = self._make_job()
        self._download()
        zip_path = os.path.join(job_dir, "result.zip")
        before = os.path.getmtime(zip_path)
        os.utime(zip_path, (before - 100, before - 100))
        self._download()
        self.assertEqual(os.path.getmtime(zip_path), before - 100)


if __name__ == "__main__":
    unittest.main()
