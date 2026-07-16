"""Tests for restrict_to_pairs on the mir-network workflow.

A network job scans the whole reference for both pools and only applies the
ceRNA pairs when the graph is built. For a pairs run that means computing the
entire genome and keeping a sliver: on job 9b4b222b, 17 of 90,515 3'UTRs and 4
of 195,159 lncRNAs -- over 99.98% of a 7h15m run discarded. restrict_to_pairs
narrows each pool to that pair side's own targets up front instead (measured:
identical predictions for every pair target, ~37s vs 7h15m).

It is opt-in because it also drops the genome-wide byproduct that /result and
the result zip would otherwise carry.

Needs the backend app deps (flask, celery, ...), so run inside the backend
container (or any env with the app_v1 stack installed):
    python3 -m unittest v2.tests.test_network_restrict_to_pairs
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


_PAIRS = [
    {"gene": "ZNRF3", "gene_refseqs": ["NM_001206998", "NM_032173"],
     "lncrna": "ENST00000842230", "score": 0.9},
    {"gene": "MTMR3", "gene_refseqs": ["NM_021090", "NM_153050"],
     "lncrna": "ENST00000602963", "score": 0.9},
    # Same gene as pair 1 -> RefSeqs must dedupe, not appear twice.
    {"gene": "ZNRF3", "gene_refseqs": ["NM_032173"],
     "lncrna": "ENST00000624012", "score": 0.9},
]


class WritePoolTargetFilesTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="isotar_restrict_")
        self._orig = app_module.BASE_DIR
        app_module.BASE_DIR = self.tmp

    def tearDown(self):
        app_module.BASE_DIR = self._orig
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _read(self, path):
        with open(path) as f:
            return [l.strip() for l in f if l.strip()]

    def test_writes_one_deduped_file_per_pool(self):
        res = app_module._write_pool_target_files("j1", _PAIRS)
        gene = self._read(res["paths"]["gene"])
        lncrna = self._read(res["paths"]["lncrna"])

        # Union across pairs, deduped (NM_032173 appears in two pairs).
        self.assertEqual(gene, ["NM_001206998", "NM_021090", "NM_032173", "NM_153050"])
        self.assertEqual(lncrna, ["ENST00000602963", "ENST00000624012", "ENST00000842230"])
        self.assertEqual(res["counts"], {"gene": 4, "lncrna": 3})

    def test_named_targets_txt_inside_each_pool_dir(self):
        # The name and location matter: -tf reads it, and result_db._build_db
        # discovers it here to filter TargetScan (which ignores the target FASTA).
        res = app_module._write_pool_target_files("j2", _PAIRS)
        output = os.path.join(app_module._job_path("j2"), "output")
        self.assertEqual(res["paths"]["gene"], os.path.join(output, "gene", "targets.txt"))
        self.assertEqual(res["paths"]["lncrna"], os.path.join(output, "lncrna", "targets.txt"))

    def test_pools_do_not_share_identifiers(self):
        # RefSeq vs Ensembl transcript -- a shared file would filter the wrong pool.
        res = app_module._write_pool_target_files("j3", _PAIRS)
        gene = self._read(res["paths"]["gene"])
        lncrna = self._read(res["paths"]["lncrna"])
        self.assertTrue(all(g.startswith("NM_") for g in gene))
        self.assertTrue(all(l.startswith("ENST") for l in lncrna))
        self.assertFalse(set(gene) & set(lncrna))


class RestrictFlagValidationTest(unittest.TestCase):
    """The flag is parsed and validated before any job state is created."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="isotar_restrict_api_")
        self._orig_base = app_module.BASE_DIR
        self._orig_default = app_module.NETWORK_RESTRICT_TO_PAIRS_DEFAULT
        self._orig_resolve = app_module._resolve_targets
        self._orig_delay = app_module.run_job.delay
        app_module.BASE_DIR = self.tmp
        app_module._resolve_targets = lambda genes, genome: ({"NM_000001"}, [])

        class _Task(object):
            id = "task-1"

        app_module.run_job.delay = lambda job_id: _Task()
        app_module.app.config["TESTING"] = True
        self.client = app_module.app.test_client()

    def tearDown(self):
        app_module.BASE_DIR = self._orig_base
        app_module.NETWORK_RESTRICT_TO_PAIRS_DEFAULT = self._orig_default
        app_module._resolve_targets = self._orig_resolve
        app_module.run_job.delay = self._orig_delay
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _post(self, payload):
        body = {
            "workflow": "mir-network",
            "mirna_ids": ["hsa-let-7a-5p"],
            "tools": ["miRanda"],
        }
        body.update(payload)
        return self.client.post("/api/v1/jobs", json=body)

    def _meta(self, resp):
        job_id = json.loads(resp.data.decode("utf-8"))["job_id"]
        with open(app_module._job_meta_path(job_id)) as f:
            return json.load(f)

    _PAIR_IN = [{"gene": "ZNRF3", "lncrna": "ENST00000842230", "score": 0.9}]

    def test_defaults_to_off_so_existing_behaviour_is_unchanged(self):
        app_module.NETWORK_RESTRICT_TO_PAIRS_DEFAULT = False
        resp = self._post({"pairs": self._PAIR_IN})
        self.assertEqual(resp.status_code, 202)
        meta = self._meta(resp)
        self.assertFalse(meta["restrict_to_pairs"])
        # No pool target files -> the runner passes None -> full reference scan.
        self.assertIsNone(meta["pool_target_files"])

    def test_opt_in_writes_pool_target_files(self):
        resp = self._post({"pairs": self._PAIR_IN, "restrict_to_pairs": True})
        self.assertEqual(resp.status_code, 202)
        meta = self._meta(resp)
        self.assertTrue(meta["restrict_to_pairs"])
        self.assertIn("gene", meta["pool_target_files"])
        self.assertIn("lncrna", meta["pool_target_files"])
        for path in meta["pool_target_files"].values():
            self.assertTrue(os.path.exists(path))

    def test_env_default_can_turn_it_on_per_deployment(self):
        app_module.NETWORK_RESTRICT_TO_PAIRS_DEFAULT = True
        resp = self._post({"pairs": self._PAIR_IN})
        self.assertEqual(resp.status_code, 202)
        self.assertTrue(self._meta(resp)["restrict_to_pairs"])

    def test_request_overrides_the_env_default(self):
        app_module.NETWORK_RESTRICT_TO_PAIRS_DEFAULT = True
        resp = self._post({"pairs": self._PAIR_IN, "restrict_to_pairs": False})
        self.assertEqual(resp.status_code, 202)
        meta = self._meta(resp)
        self.assertFalse(meta["restrict_to_pairs"])
        self.assertIsNone(meta["pool_target_files"])

    def test_restrict_without_pairs_is_rejected(self):
        # Nothing to restrict to; silently running genome-wide would be a lie.
        resp = self._post({"restrict_to_pairs": True})
        self.assertEqual(resp.status_code, 400)
        self.assertIn("requires pairs", json.loads(resp.data.decode("utf-8"))["error"])

    def test_non_boolean_is_rejected(self):
        resp = self._post({"pairs": self._PAIR_IN, "restrict_to_pairs": "yes"})
        self.assertEqual(resp.status_code, 400)

    def test_discovery_run_unaffected(self):
        # No pairs, no flag -> full reference, exactly as before.
        resp = self._post({})
        self.assertEqual(resp.status_code, 202)
        meta = self._meta(resp)
        self.assertFalse(meta["restrict_to_pairs"])
        self.assertIsNone(meta["pool_target_files"])


class BuildDbTargetsDiscoveryTest(unittest.TestCase):
    """result_db._build_db must find a pool-local targets.txt (restrict_to_pairs
    writes one per pool) while single-pool jobs keep using the job-dir one."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="isotar_targets_discovery_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _resolve(self, output_dir):
        """Mirror the lookup in _build_db."""
        candidate = os.path.join(output_dir, "targets.txt")
        if not os.path.exists(candidate):
            candidate = os.path.join(
                os.path.dirname(os.path.abspath(output_dir)), "targets.txt")
        return candidate

    def test_pool_local_file_wins(self):
        pool_dir = os.path.join(self.tmp, "output", "gene")
        os.makedirs(pool_dir)
        pool_file = os.path.join(pool_dir, "targets.txt")
        open(pool_file, "w").write("NM_1\n")
        # A stale job-dir file must not shadow the pool's own list.
        open(os.path.join(self.tmp, "output", "targets.txt"), "w").write("NM_999\n")
        self.assertEqual(self._resolve(pool_dir), pool_file)

    def test_falls_back_to_job_dir_for_single_pool_jobs(self):
        output_dir = os.path.join(self.tmp, "output")
        os.makedirs(output_dir)
        job_file = os.path.join(self.tmp, "targets.txt")
        open(job_file, "w").write("NM_1\n")
        self.assertEqual(self._resolve(output_dir), job_file)


if __name__ == "__main__":
    unittest.main()
