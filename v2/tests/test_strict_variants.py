"""mirna_processing.py --strict: a requested variant that cannot be built is an
error, not a skipped warning.

Job b6ddbc3a asked for 4:C|U on hsa-miR-101-3p (UACAGUACUGUGAUAACUGAA), whose
position 4 is A: without --strict the script skipped the modification and
wrote the wild type alone, and the job reported success.

Runs the script end to end, as app_v1 does, under the current interpreter and
under python2.7 when it is installed (the script must stay 2.7-compatible).

Run with:
    python3 -m unittest v2.tests.test_strict_variants
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SCRIPT = os.path.join(_REPO_ROOT, "v2", "mirna_processing.py")
SHIPPED_RESOURCES_DIR = os.path.join(_REPO_ROOT, "v2", "opt", "resources")

MIRNA = "hsa-miR-101-3p"
PRE = "hsa-mir-101-2"
MATURE = "UACAGUACUGUGAUAACUGAA"


def _interpreters():
    found = [sys.executable]
    py27 = shutil.which("python2.7")
    if py27:
        found.append(py27)
    return found


@unittest.skipUnless(
    os.path.exists(os.path.join(SHIPPED_RESOURCES_DIR, "mature_pre_mirna_ext.json")),
    "shipped resources dir not present (expected at v2/opt/resources)")
class StrictVariantTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.env = dict(os.environ, ISOTAR_RESOURCES_DIR=SHIPPED_RESOURCES_DIR)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, python, *extra):
        out = os.path.join(self.tmp, "mirna.fa")
        if os.path.exists(out):
            os.remove(out)
        proc = subprocess.run(
            [python, SCRIPT, MIRNA, "-o", out, "--pre-id", PRE] + list(extra),
            env=self.env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True)
        headers = []
        if os.path.exists(out):
            with open(out) as f:
                headers = [l[1:].strip() for l in f if l.startswith(">")]
        return proc, headers

    def test_mismatched_base_is_an_error_with_strict(self):
        for python in _interpreters():
            with self.subTest(python=python):
                proc, headers = self._run(python, "-m", "4:C|U", "--strict")
                self.assertEqual(proc.returncode, 1)
                self.assertIn("Expected 'C' at position 4, found 'A'", proc.stderr)
                self.assertIn("4:C|U", proc.stderr)
                self.assertEqual(headers, [])  # nothing written

    def test_mismatched_base_is_still_skipped_without_strict(self):
        # Legacy behaviour kept for direct CLI use.
        proc, headers = self._run(sys.executable, "-m", "4:C|U")
        self.assertEqual(proc.returncode, 0)
        self.assertEqual(headers, ["{},WT".format(MIRNA)])
        self.assertIn("Skipping modification", proc.stdout)

    def test_matching_base_builds_the_variant_with_strict(self):
        for python in _interpreters():
            with self.subTest(python=python):
                proc, headers = self._run(python, "-m", "3:C|U", "--strict")
                self.assertEqual(proc.returncode, 0, proc.stderr)
                self.assertEqual(headers, ["{},WT".format(MIRNA),
                                           "{},3:C|U,modified".format(MIRNA)])

    def test_variant_out_of_length_range_is_an_error_with_strict(self):
        # -10|0 extends the 21 nt mature to 31 nt, past the 30 nt limit.
        proc, headers = self._run(sys.executable, "--shift=-10|0", "--strict")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("31 nt", proc.stderr)
        proc, headers = self._run(sys.executable, "--shift=-10|0")
        self.assertEqual(proc.returncode, 0)
        self.assertEqual(headers, ["{},WT".format(MIRNA)])  # dropped, legacy

    def test_combined_mismatch_is_an_error_with_strict(self):
        # After a 0|0 shift the sequence is the mature itself: position 4 is A.
        proc, _ = self._run(sys.executable, "-m", "4:C|U", "--shift=0|0", "-b", "--strict")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Expected 'C' at position 4, found 'A'", proc.stderr)


if __name__ == "__main__":
    unittest.main()
