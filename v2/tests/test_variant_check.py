"""Tests for app_v1.variant_check: requested variants missing from mirna.fa.

Run with:
    python3 -m unittest v2.tests.test_variant_check
"""

import os
import sys
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app_v1.variant_check import (  # noqa: E402
    applied_variants, requested_variants, variant_warnings)

MIR = "hsa-miR-101-3p"


def single(mods=None, shift=None, **extra):
    meta = {"workflow": "mir-target", "mirna_id": MIR, "mirna_seq": None,
            "modifications": mods or [], "shift": shift}
    meta.update(extra)
    return meta


class VariantCheckTests(unittest.TestCase):
    def test_job_b6ddbc3a_ran_wild_type_only(self):
        warnings = variant_warnings(single(["4:C|U"]), [">{},WT\n".format(MIR)])
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(w["code"], "variant_not_applied")
        self.assertEqual(w["mirna_id"], MIR)
        self.assertEqual(w["missing"], {"modifications": ["4:C|U"], "shift": None})
        self.assertIn("wild type only", w["message"])
        self.assertIn("4:C|U", w["message"])

    def test_applied_modification_gives_no_warning(self):
        headers = [">{},WT".format(MIR), ">{},3:C|U,modified".format(MIR)]
        self.assertEqual(variant_warnings(single(["3:C|U"]), headers), [])

    def test_partly_applied_modifications(self):
        headers = [">{},WT".format(MIR), ">{},3:C|U,modified".format(MIR)]
        warnings = variant_warnings(single(["3:C|U", "4:C|U"]), headers)
        self.assertEqual(warnings[0]["missing"]["modifications"], ["4:C|U"])
        self.assertIn("Part of the requested variant", warnings[0]["message"])

    def test_shift_is_compared_numerically(self):
        headers = [">{},WT".format(MIR), ">{},1|0,shifted".format(MIR)]
        self.assertEqual(variant_warnings(single(shift="+1|0"), headers), [])
        self.assertEqual(len(variant_warnings(single(shift="-1|0"), headers)), 1)

    def test_combined_operation(self):
        headers = [">{},WT".format(MIR), ">{},2:G|A&5:A|U,-7|1,modified_shifted".format(MIR)]
        meta = single(["2:G|A", "5:A|U"], "-7|1")
        self.assertEqual(variant_warnings(meta, headers), [])
        applied = applied_variants(headers)[MIR]
        self.assertEqual(applied["modifications"], {"2:G|A", "5:A|U"})
        self.assertEqual(applied["shifts"], {(-7, 1)})

    def test_no_request_and_custom_sequence_give_no_warning(self):
        self.assertEqual(requested_variants(single()), {})
        self.assertEqual(requested_variants(single(["4:C|U"], mirna_seq="UACAGU")), {})
        self.assertEqual(variant_warnings(single(), [">{},WT".format(MIR)]), [])

    def test_network_variants(self):
        meta = {"workflow": "mir-network", "mirna_ids": ["hsa-miR-376a-5p"],
                "mirna_variants": {"hsa-miR-376a-5p": [
                    {"shift": None, "modifications": ["4:A|G"]},
                    {"shift": "1|0", "modifications": []}]}}
        headers = [">hsa-miR-376a-5p,WT", ">hsa-miR-376a-5p,1|0,shifted"]
        warnings = variant_warnings(meta, headers)
        self.assertEqual(len(warnings), 1)  # the shift ran, the modification did not
        self.assertEqual(warnings[0]["missing"]["modifications"], ["4:A|G"])


if __name__ == "__main__":
    unittest.main()
