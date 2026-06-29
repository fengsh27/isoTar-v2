"""Tests for per-miRNA precursor selection in the miR-Network workflow.

Covers the two halves of the `pre_ids` feature in app_v1.app:
  1. Request validation in _submit_network_job (the POST /api/v1/jobs handler)
     -- pre_ids must be an object of non-empty id->pre_id strings whose keys all
     reference a submitted miRNA, and a valid map is persisted into job.json.
  2. Command building in _process_mirna_list -- a mapped miRNA gets `--pre-id`
     appended to its mirna_processing.py invocation; an unmapped one does not.

Needs the backend app deps (flask, celery, ...), so run inside the backend
container (or any env with the app_v1 stack installed):
    python3 -m unittest v2.tests.test_network_pre_ids
"""

import json
import os
import shutil
import sys
import tempfile
import types
import unittest
from unittest import mock

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app_v1 import app as app_module  # noqa: E402


class _BaseDirTestCase(unittest.TestCase):
    """Redirect job state into a temp dir so tests never touch /opt/out/jobs."""

    def setUp(self):
        self._tmp = tempfile.mkdtemp(prefix="isotar_test_jobs_")
        self._orig_base_dir = app_module.BASE_DIR
        # _job_path() reads BASE_DIR at call time, so overriding the module
        # global is enough to reroute every job/meta path below it.
        app_module.BASE_DIR = self._tmp

    def tearDown(self):
        app_module.BASE_DIR = self._orig_base_dir
        shutil.rmtree(self._tmp, ignore_errors=True)


class PreIdsValidationTests(_BaseDirTestCase):
    """Validation of the optional `pre_ids` field on a mir-network submission."""

    def setUp(self):
        super(PreIdsValidationTests, self).setUp()
        self.client = app_module.app.test_client()

    def _post(self, body):
        return self.client.post("/api/v1/jobs", json=body)

    def test_rejects_non_dict_pre_ids(self):
        resp = self._post({
            "workflow": "mir-network",
            "mirna_ids": ["hsa-let-7a-5p"],
            "tools": ["miRanda"],
            "pre_ids": "hsa-let-7a-1",  # not an object
        })
        self.assertEqual(resp.status_code, 400)
        self.assertIn("pre_ids must be an object", resp.get_json()["error"])

    def test_rejects_empty_pre_id_value(self):
        resp = self._post({
            "workflow": "mir-network",
            "mirna_ids": ["hsa-let-7a-5p"],
            "tools": ["miRanda"],
            "pre_ids": {"hsa-let-7a-5p": "  "},  # blank precursor id
        })
        self.assertEqual(resp.status_code, 400)
        self.assertIn("non-empty", resp.get_json()["error"])

    def test_rejects_pre_id_key_not_in_mirna_ids(self):
        resp = self._post({
            "workflow": "mir-network",
            "mirna_ids": ["hsa-miR-21-5p"],
            "tools": ["miRanda"],
            "pre_ids": {"hsa-let-7a-5p": "hsa-let-7a-1"},  # stray key (typo guard)
        })
        self.assertEqual(resp.status_code, 400)
        body = resp.get_json()
        self.assertIn("not in mirna_ids", body["error"])
        self.assertEqual(body["unknown"], ["hsa-let-7a-5p"])

    def test_valid_pre_ids_persisted_in_meta(self):
        # Mock the Celery dispatch so no broker is needed and run_job never fires.
        fake_task = types.SimpleNamespace(id="task-test-123")
        with mock.patch.object(app_module.run_job, "delay", return_value=fake_task) as delay:
            resp = self._post({
                "workflow": "mir-network",
                "mirna_ids": ["hsa-let-7a-5p", "hsa-miR-21-5p"],
                "tools": ["miRanda"],
                "cores": 1,
                "pre_ids": {"hsa-let-7a-5p": "hsa-let-7a-1"},
            })
        self.assertEqual(resp.status_code, 202)
        delay.assert_called_once()
        job_id = resp.get_json()["job_id"]

        with open(os.path.join(self._tmp, job_id, "job.json")) as f:
            meta = json.load(f)
        self.assertEqual(meta["pre_ids"], {"hsa-let-7a-5p": "hsa-let-7a-1"})

    def test_pre_ids_omitted_defaults_to_empty(self):
        # Back-compat: a request with no pre_ids still succeeds and stores {}.
        fake_task = types.SimpleNamespace(id="task-test-456")
        with mock.patch.object(app_module.run_job, "delay", return_value=fake_task):
            resp = self._post({
                "workflow": "mir-network",
                "mirna_ids": ["hsa-miR-21-5p"],
                "tools": ["miRanda"],
                "cores": 1,
            })
        self.assertEqual(resp.status_code, 202)
        job_id = resp.get_json()["job_id"]
        with open(os.path.join(self._tmp, job_id, "job.json")) as f:
            meta = json.load(f)
        self.assertEqual(meta["pre_ids"], {})


class ProcessMirnaListPreIdTests(_BaseDirTestCase):
    """_process_mirna_list passes --pre-id only for miRNAs present in pre_ids."""

    def test_pre_id_appended_only_for_mapped_mirna(self):
        job_id = "job-preid"
        os.makedirs(app_module._job_path(job_id))
        mirna_ids = ["hsa-let-7a-5p", "hsa-miR-21-5p"]
        meta = {"pre_ids": {"hsa-let-7a-5p": "hsa-let-7a-1"}}

        captured = []

        def fake_check_call(cmd, **kwargs):
            captured.append(cmd)
            # Write a stub FASTA to the -o path so the concat step has input.
            out_path = cmd[cmd.index("-o") + 1]
            with open(out_path, "w") as f:
                f.write(">{},WT\nACGUACGUACGUACGUACGU\n".format(cmd[2]))
            return 0

        fasta_path = os.path.join(app_module._job_path(job_id), "mirna.fa")
        with mock.patch.object(app_module.subprocess, "check_call", side_effect=fake_check_call):
            app_module._process_mirna_list(job_id, mirna_ids, fasta_path, meta)

        cmd_by_mirna = {cmd[2]: cmd for cmd in captured}
        self.assertEqual(set(cmd_by_mirna), set(mirna_ids))

        # Mapped miRNA -> '--pre-id hsa-let-7a-1' present, as adjacent tokens.
        mapped = cmd_by_mirna["hsa-let-7a-5p"]
        self.assertIn("--pre-id", mapped)
        self.assertEqual(mapped[mapped.index("--pre-id") + 1], "hsa-let-7a-1")

        # Unmapped miRNA -> no precursor flag.
        self.assertNotIn("--pre-id", cmd_by_mirna["hsa-miR-21-5p"])

        # No spurious errors recorded; concatenated FASTA written.
        self.assertNotIn("mirna_errors", meta)
        self.assertTrue(os.path.exists(fasta_path))


if __name__ == "__main__":
    unittest.main()
