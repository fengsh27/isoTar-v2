"""Tests for TargetScan filename sanitization (v2.mirna_predicting).

TargetScan's vendored perl (targetscan_70_BL_PCT.pl) shells out
`sort -k8,8n FILE > FILE.sort.txt` with the filename *unquoted*. A shift/
modification variant header carries '|' (shift/mod separator) and '&' (combined-
mod joiner), both shell-special -- so feeding the raw header into TargetScan's
working filenames breaks the sort (perl dies, the tool is left half-run). We
give the perl a sanitized basename while keeping the raw header on the final
merged results (Python writes those; parse_result.py reads them back by header).

The unit tests run anywhere. The end-to-end test needs the TargetScan install +
datasets, so it runs in the backend container and skips elsewhere:
    python3 -m unittest v2.tests.test_targetscan_safe_names
"""

import os
import shutil
import sys
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from v2.mirna_predicting import (  # noqa: E402
    _safe_ts_basename,
    targetscan_prep,
    run_targetscan,
)

# TargetScan needs its own precomputed 3' UTR + conservation datasets; present
# only where TargetScan is installed (the backend image).
_TS_UTR0 = "/opt/TargetScan/Datasets/3utr/targetscan_utr_part_0.txt"
_TS_BLN0 = "/opt/TargetScan/Datasets/bln_bins/targetscan_median_bls_bins_part_0.txt"
_TS_FAMILY = "/opt/TargetScan/Datasets/miR_Family_Info.json"
_HAVE_TS = all(os.path.exists(p) for p in (_TS_UTR0, _TS_BLN0, _TS_FAMILY))

# The exact variant header that crashed a live job (see investigation).
_VARIANT_HEADER = "hsa-let-7a-2-3p,2|3,modified_shifted"
_VARIANT_SEQ = "GUACAGCCUCCUAGCUUUCCUUG"


class SafeBasenameTests(unittest.TestCase):
    def test_pipe_and_amp_and_comma_replaced(self):
        self.assertEqual(
            _safe_ts_basename("hsa-let-7a-2-3p,2|3,modified_shifted"),
            "hsa-let-7a-2-3p_2_3_modified_shifted",
        )

    def test_combined_mod_ampersand_replaced(self):
        # Combined modifications join with '&' (also shell-special).
        self.assertEqual(
            _safe_ts_basename("hsa-miR-21-5p,8:A|U&2:C|G,modified"),
            "hsa-miR-21-5p_8_A_U_2_C_G_modified",
        )

    def test_safe_chars_preserved(self):
        # Letters, digits, dot, underscore, hyphen must survive unchanged.
        self.assertEqual(_safe_ts_basename("hsa-let-7a-2-3p,WT"),
                         "hsa-let-7a-2-3p_WT")

    def test_no_shell_metacharacters_remain(self):
        for ch in "|&,;<>()$`\\\"' \t*?":
            self.assertNotIn(ch, _safe_ts_basename("x{}y".format(ch)))

    def test_wt_and_variant_stay_distinct(self):
        # Sanitization must not collapse a miRNA's WT and variant to one name,
        # or their TargetScan results (distinguished only by filename) would merge.
        self.assertNotEqual(
            _safe_ts_basename("hsa-let-7a-2-3p,WT"),
            _safe_ts_basename("hsa-let-7a-2-3p,2|3,modified_shifted"),
        )


@unittest.skipUnless(_HAVE_TS, "TargetScan install/datasets not present")
class TargetscanVariantHeaderE2ETests(unittest.TestCase):
    """The real perl must run to completion on a '|'-bearing variant header, and
    must NOT emit the shell-pipe artifact (a file named after the header
    truncated at the first '|')."""

    def setUp(self):
        self._dir = tempfile.mkdtemp(prefix="ts_safe_")

    def tearDown(self):
        shutil.rmtree(self._dir, ignore_errors=True)

    def test_variant_header_runs_and_leaves_no_shell_artifact(self):
        safe = _safe_ts_basename(_VARIANT_HEADER)
        targetscan_prep(_VARIANT_SEQ, _VARIANT_HEADER, self._dir)
        # prep must name the input by the SAFE basename, not the raw header.
        ts_input = os.path.join(self._dir, "{}_targetscan.txt".format(safe))
        self.assertTrue(os.path.exists(ts_input),
                        "targetscan_prep should write a shell-safe input filename")

        out1 = os.path.join(self._dir, "{}_part_0_out1.txt".format(safe))
        out2 = os.path.join(self._dir, "{}_part_0_out2.txt".format(safe))
        # Before the fix this raised CalledProcessError (perl died on the pipe).
        run_targetscan(ts_input, _TS_UTR0, out1, _TS_BLN0, out2)

        self.assertTrue(os.path.exists(out1) and os.path.getsize(out1) > 0)
        self.assertTrue(os.path.exists(out2))

        # The shell-pipe fingerprint: a file named after the header truncated at
        # the first '|'. It must NOT exist.
        artifact = os.path.join(self._dir, _VARIANT_HEADER.split("|")[0])
        self.assertFalse(os.path.exists(artifact),
                         "shell-pipe artifact leaked: {}".format(artifact))

        # The final merged results legitimately keep the raw '|' header -- Python
        # writes them, so no shell is involved. Confirm that still works.
        merged = os.path.join(self._dir,
                              "{}_Targetscan_results1.txt".format(_VARIANT_HEADER))
        with open(out1) as src, open(merged, "w") as dst:
            dst.write(src.read())
        self.assertTrue(os.path.exists(merged))


if __name__ == "__main__":
    unittest.main()
