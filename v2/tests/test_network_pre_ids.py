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

    def test_valid_shifts_and_modifications_persisted(self):
        fake_task = types.SimpleNamespace(id="task-ops")
        with mock.patch.object(app_module.run_job, "delay", return_value=fake_task):
            resp = self._post({
                "workflow": "mir-network",
                "mirna_ids": ["hsa-miR-21-5p", "hsa-miR-155-5p"],
                "tools": ["miRanda"],
                "cores": 1,
                "shifts": {"hsa-miR-21-5p": "-7|1"},
                "modifications": {"hsa-miR-155-5p": ["8:A>U"]},
            })
        self.assertEqual(resp.status_code, 202)
        job_id = resp.get_json()["job_id"]
        with open(os.path.join(self._tmp, job_id, "job.json")) as f:
            meta = json.load(f)
        self.assertEqual(meta["shifts"], {"hsa-miR-21-5p": "-7|1"})
        self.assertEqual(meta["modifications"], {"hsa-miR-155-5p": ["8:A>U"]})

    def test_shifts_referencing_unknown_mirna_rejected(self):
        resp = self._post({
            "workflow": "mir-network",
            "mirna_ids": ["hsa-miR-21-5p"],
            "tools": ["miRanda"],
            "shifts": {"hsa-miR-999-5p": "-7|1"},
        })
        self.assertEqual(resp.status_code, 400)
        self.assertIn("hsa-miR-999-5p", resp.get_json()["unknown"])

    def test_modifications_must_be_lists(self):
        resp = self._post({
            "workflow": "mir-network",
            "mirna_ids": ["hsa-miR-21-5p"],
            "tools": ["miRanda"],
            "modifications": {"hsa-miR-21-5p": "8:A>U"},  # string, not a list
        })
        self.assertEqual(resp.status_code, 400)

    def test_variants_normalized_and_persisted(self):
        fake_task = types.SimpleNamespace(id="task-var")
        with mock.patch.object(app_module.run_job, "delay", return_value=fake_task):
            resp = self._post({
                "workflow": "mir-network",
                "mirna_ids": ["hsa-miR-21-5p"],
                "tools": ["miRanda"],
                "cores": 1,
                "variants": {"hsa-miR-21-5p": [
                    {"shift": "-7|1"}, {"shift": "-4|0"}, {"modifications": ["8:A|U"]},
                ]},
            })
        self.assertEqual(resp.status_code, 202)
        job_id = resp.get_json()["job_id"]
        with open(os.path.join(self._tmp, job_id, "job.json")) as f:
            meta = json.load(f)
        specs = meta["mirna_variants"]["hsa-miR-21-5p"]
        self.assertEqual(len(specs), 3)
        self.assertEqual(specs[0], {"shift": "-7|1", "modifications": []})
        self.assertEqual(specs[2], {"shift": None, "modifications": ["8:A|U"]})

    def test_variants_merge_with_single_maps(self):
        # shifts/modifications fold in as an implicit spec alongside variants.
        fake_task = types.SimpleNamespace(id="task-var2")
        with mock.patch.object(app_module.run_job, "delay", return_value=fake_task):
            resp = self._post({
                "workflow": "mir-network",
                "mirna_ids": ["hsa-miR-21-5p"],
                "tools": ["miRanda"],
                "cores": 1,
                "shifts": {"hsa-miR-21-5p": "-7|1"},
                "variants": {"hsa-miR-21-5p": [{"modifications": ["8:A|U"]}]},
            })
        self.assertEqual(resp.status_code, 202)
        job_id = resp.get_json()["job_id"]
        with open(os.path.join(self._tmp, job_id, "job.json")) as f:
            meta = json.load(f)
        specs = meta["mirna_variants"]["hsa-miR-21-5p"]
        self.assertEqual(len(specs), 2)
        self.assertIn({"shift": "-7|1", "modifications": []}, specs)
        self.assertIn({"shift": None, "modifications": ["8:A|U"]}, specs)

    def test_variant_referencing_unknown_mirna_rejected(self):
        resp = self._post({
            "workflow": "mir-network",
            "mirna_ids": ["hsa-miR-21-5p"],
            "tools": ["miRanda"],
            "variants": {"hsa-miR-999-5p": [{"shift": "-7|1"}]},
        })
        self.assertEqual(resp.status_code, 400)
        self.assertIn("hsa-miR-999-5p", resp.get_json()["unknown"])

    def test_empty_variant_spec_rejected(self):
        resp = self._post({
            "workflow": "mir-network",
            "mirna_ids": ["hsa-miR-21-5p"],
            "tools": ["miRanda"],
            "variants": {"hsa-miR-21-5p": [{}]},  # neither shift nor modifications
        })
        self.assertEqual(resp.status_code, 400)


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

    def test_shift_and_modifications_appended_per_mirna(self):
        job_id = "job-ops"
        os.makedirs(app_module._job_path(job_id))
        mirna_ids = ["hsa-miR-21-5p", "hsa-let-7a-5p"]
        meta = {
            "shifts": {"hsa-miR-21-5p": "-7|1"},
            "modifications": {"hsa-miR-21-5p": ["8:A>U"]},
        }

        captured = []

        def fake_check_call(cmd, **kwargs):
            captured.append(cmd)
            out_path = cmd[cmd.index("-o") + 1]
            with open(out_path, "w") as f:
                f.write(">{},WT\nACGUACGUACGUACGUACGU\n".format(cmd[2]))
            return 0

        fasta_path = os.path.join(app_module._job_path(job_id), "mirna.fa")
        with mock.patch.object(app_module.subprocess, "check_call", side_effect=fake_check_call):
            app_module._process_mirna_list(job_id, mirna_ids, fasta_path, meta)

        by_mirna = {cmd[2]: cmd for cmd in captured}
        m21 = by_mirna["hsa-miR-21-5p"]
        self.assertIn("-m", m21)
        self.assertEqual(m21[m21.index("-m") + 1], "8:A>U")
        self.assertIn("--shift=-7|1", m21)   # =-form for the negative shift
        self.assertIn("-b", m21)             # both mods and shift -> combined

        # The miRNA with no ops gets none of the operation flags.
        let7 = by_mirna["hsa-let-7a-5p"]
        self.assertNotIn("-m", let7)
        self.assertNotIn("-b", let7)
        self.assertFalse(any(a.startswith("--shift=") for a in let7))

    def test_multiple_variants_run_once_each_and_dedup_wt(self):
        job_id = "job-multivar"
        os.makedirs(app_module._job_path(job_id))
        mirna_ids = ["hsa-miR-21-5p"]
        meta = {"mirna_variants": {"hsa-miR-21-5p": [
            {"shift": "-7|1", "modifications": []},
            {"shift": "-4|0", "modifications": []},
        ]}}

        captured = []

        def fake_check_call(cmd, **kwargs):
            captured.append(cmd)
            out_path = cmd[cmd.index("-o") + 1]
            mid = cmd[2]
            # Emulate mirna_processing: always a WT record, plus one variant.
            recs = [">{},WT".format(mid), "ACGUACGUACGUACGUACGU"]
            if "--shift=-7|1" in cmd:
                recs += [">{},-7|1,shifted".format(mid), "ACGUACGUACGUACGUAC"]
            elif "--shift=-4|0" in cmd:
                recs += [">{},-4|0,shifted".format(mid), "ACGUACGUACGUACGUACGUAC"]
            with open(out_path, "w") as f:
                f.write("\n".join(recs) + "\n")
            return 0

        fasta_path = os.path.join(app_module._job_path(job_id), "mirna.fa")
        with mock.patch.object(app_module.subprocess, "check_call", side_effect=fake_check_call):
            app_module._process_mirna_list(job_id, mirna_ids, fasta_path, meta)

        # One subprocess run per variant spec.
        self.assertEqual(len(captured), 2)
        # Concatenated FASTA: a single (deduped) WT + one record per variant.
        with open(fasta_path) as f:
            headers = [ln.strip() for ln in f if ln.startswith(">")]
        self.assertEqual(headers.count(">hsa-miR-21-5p,WT"), 1)
        self.assertIn(">hsa-miR-21-5p,-7|1,shifted", headers)
        self.assertIn(">hsa-miR-21-5p,-4|0,shifted", headers)
        self.assertEqual(len(headers), 3)


if __name__ == "__main__":
    unittest.main()
