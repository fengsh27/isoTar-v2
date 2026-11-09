"""Regression tests for the mature-window coordinate convention.

The shipped catalogs do not agree on how `ext_mature_loc_start` is counted.
The un-prefixed human file stores it 0-based (its mature region is
`ext_pre_seq[start:end]`); the nine `<sp>_mature_pre_mirna_ext.json` files
store it 1-based inclusive. `ext_mature_loc_end` names the same position under
both conventions, so only the start differs, and the field name gives no hint
which a file used.

`apply_shift` used to slice 0-based. That matched the human file and nothing
else, so on all nine non-human species a shifted variant was analysed one base
3' of its true window and one nt short -- losing the 5' base that sets the seed
register (positions 2-8). The old bound checks (`0 <= start < end <= len(pre)`)
pass under either convention, which is why it went unnoticed; the tests here
tie the numbers to the letters instead.

`find_mirna_sequence` now returns a 1-based inclusive start, normalized by
`anchored_mature_start`, and `apply_shift` slices to match.

Runnable with:  python3 -m unittest v2.tests.test_mature_coordinates
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from v2.mirna_processing import (  # noqa: E402
    LEGACY_HUMAN_FILE,
    anchored_mature_start,
    apply_shift,
    find_mirna_sequence,
    load_json,
    resolve_metadata_path,
)

SHIPPED_RESOURCES_DIR = os.path.join(_REPO_ROOT, "v2", "opt", "resources")

SHIPPED_NON_HUMAN_SPECIES = (
    "cel", "cfa", "dme", "dre", "mdo", "mml", "mmu", "ptr", "rno",
)

# Known-good reference points, straight from miRBase. Both files store the same
# mature sequence for this miRNA, under opposite conventions: the human record
# is ext_mature_loc_start=157 (0-based), the mouse one 68 (1-based).
MIR21 = "UAGCUUAUCAGACUGAUGUUGA"


def _first_clean_id(json_path):
    """First id in `json_path` with a single, fully-populated precursor.

    Skips the handful of entries with an empty `ext_pre_seq` and the ids that
    map to several precursors and so need --pre-id, neither of which is what
    these tests are about.
    """
    with open(json_path, "r") as f:
        data = json.load(f)
    for mirna_id, entries in sorted(data.items()):
        if not entries or len(entries) != 1:
            continue
        e = entries[0]
        if (e.get("mature_seq") and e.get("ext_pre_seq")
                and isinstance(e.get("ext_mature_loc_start"), int)
                and isinstance(e.get("ext_mature_loc_end"), int)):
            return mirna_id
    raise AssertionError(
        "no single-precursor, fully-populated id in {}".format(json_path)
    )


def _shipped_files():
    return sorted(glob.glob(
        os.path.join(SHIPPED_RESOURCES_DIR, "*mature_pre_mirna_ext.json")
    ))


class AnchoredMatureStartTest(unittest.TestCase):
    """The convention is decided per record, by comparison with mature_seq."""

    def test_keeps_a_one_based_start(self):
        pre = "GGG" + MIR21 + "CCC"          # mature at 0-based index 3
        self.assertEqual(anchored_mature_start(pre, MIR21, 4), 4)

    def test_promotes_a_zero_based_start(self):
        pre = "GGG" + MIR21 + "CCC"
        self.assertEqual(anchored_mature_start(pre, MIR21, 3), 4)

    def test_prefers_the_stored_position_when_the_mature_repeats(self):
        # 19 shipped records contain mature_seq more than once (hsa-miR-941
        # five times), so a blind search could pick the wrong copy. The stored
        # value is tested first precisely so it wins here.
        pre = MIR21 + "AAAA" + MIR21
        second = len(MIR21) + 4 + 1          # 1-based start of the 2nd copy
        self.assertEqual(anchored_mature_start(pre, MIR21, second), second)
        self.assertEqual(anchored_mature_start(pre, MIR21, second - 1), second)

    def test_searches_when_the_stored_start_is_unusable(self):
        pre = "GGG" + MIR21 + "CCC"
        self.assertEqual(anchored_mature_start(pre, MIR21, 999), 4)

    def test_returns_the_stored_start_when_the_mature_is_absent(self):
        self.assertEqual(anchored_mature_start("GGGGCCCC", MIR21, 7), 7)

    def test_returns_the_stored_start_for_an_empty_mature(self):
        self.assertEqual(anchored_mature_start("GGGG", "", 2), 2)


@unittest.skipUnless(
    os.path.isdir(SHIPPED_RESOURCES_DIR),
    "shipped resources dir not present (expected at v2/opt/resources)",
)
class ShippedCatalogWindowTest(unittest.TestCase):
    """Sweep every record in every shipped catalog, not just a sample."""

    def test_every_record_resolves_to_its_mature_sequence(self):
        checked, files = 0, 0
        for path in _shipped_files():
            files += 1
            with open(path) as f:
                data = json.load(f)
            for mirna_id, entries in data.items():
                for entry in entries:
                    pre = (entry.get("ext_pre_seq") or "").upper()
                    mature = (entry.get("mature_seq") or "").upper()
                    if not pre or not mature:
                        continue
                    stored = entry["ext_mature_loc_start"]
                    start = anchored_mature_start(pre, mature, stored)
                    end = entry["ext_mature_loc_end"]
                    self.assertEqual(
                        pre[start - 1:end], mature,
                        "{} in {}: window {}-{} (stored start {}) is not the "
                        "mature sequence".format(
                            mirna_id, os.path.basename(path), start, end, stored,
                        ),
                    )
                    checked += 1
        self.assertEqual(files, 10, "expected 10 shipped catalogs")
        self.assertGreater(checked, 9000,
                           "only {} records checked".format(checked))

    def test_the_stored_end_needs_no_adjustment(self):
        """A substring at 0-based [i, i+L) sits at 1-based i+1 .. i+L, so the
        exclusive 0-based end and the inclusive 1-based end are the same
        number. This is why only the start ever differed."""
        for path in _shipped_files():
            with open(path) as f:
                data = json.load(f)
            for mirna_id, entries in data.items():
                for entry in entries:
                    pre = (entry.get("ext_pre_seq") or "").upper()
                    mature = (entry.get("mature_seq") or "").upper()
                    if not pre or not mature:
                        continue
                    index = pre.find(mature)
                    if index < 0:
                        continue
                    self.assertEqual(
                        entry["ext_mature_loc_end"], index + len(mature),
                        "{} in {}: unexpected end coordinate".format(
                            mirna_id, os.path.basename(path),
                        ),
                    )

    def test_both_conventions_are_present_in_the_shipped_data(self):
        """Guards the premise: if a future regeneration unifies the files,
        this fails and the comment in anchored_mature_start can be retired."""
        human = os.path.join(SHIPPED_RESOURCES_DIR, LEGACY_HUMAN_FILE)
        with open(human) as f:
            entry = json.load(f)["hsa-miR-21-5p"][0]
        self.assertEqual(entry["mature_seq"], MIR21)
        self.assertEqual(entry["ext_mature_loc_start"], 157)   # 0-based
        self.assertEqual(entry["ext_pre_seq"].find(MIR21), 157)

        mouse = os.path.join(SHIPPED_RESOURCES_DIR,
                             "mmu_mature_pre_mirna_ext.json")
        with open(mouse) as f:
            entry = json.load(f)["mmu-miR-21a-5p"][0]
        self.assertEqual(entry["mature_seq"], MIR21)
        self.assertEqual(entry["ext_mature_loc_start"], 68)    # 1-based
        self.assertEqual(entry["ext_pre_seq"].find(MIR21), 67)


@unittest.skipUnless(
    os.path.isdir(SHIPPED_RESOURCES_DIR),
    "shipped resources dir not present (expected at v2/opt/resources)",
)
class ApplyShiftTest(unittest.TestCase):
    """What the prediction tools actually receive."""

    def _resolve(self, fname):
        mirna_id = _first_clean_id(os.path.join(SHIPPED_RESOURCES_DIR, fname))
        path = resolve_metadata_path(mirna_id,
                                     resources_dir=SHIPPED_RESOURCES_DIR)
        mature, pre, start, end, _ = find_mirna_sequence(load_json(path), mirna_id)
        return mirna_id, mature, pre, start, end

    def _all_species_files(self):
        return [LEGACY_HUMAN_FILE] + [
            "{}_mature_pre_mirna_ext.json".format(s)
            for s in SHIPPED_NON_HUMAN_SPECIES
        ]

    def test_zero_shift_is_exactly_the_mature_sequence(self):
        """The bug in one assertion: a 0|0 variant is the wild type."""
        for fname in self._all_species_files():
            with self.subTest(catalog=fname):
                mirna_id, mature, pre, start, end = self._resolve(fname)
                shifted, info = apply_shift(pre, mature, start, end, "0|0")
                self.assertEqual(shifted, mature,
                                 "0|0 shift changed the sequence for {}"
                                 .format(mirna_id))
                self.assertEqual(info, "0|0")

    def test_left_shift_trims_only_the_five_prime_end(self):
        for fname in self._all_species_files():
            with self.subTest(catalog=fname):
                _, mature, pre, start, end = self._resolve(fname)
                shifted, _ = apply_shift(pre, mature, start, end, "1|0")
                self.assertEqual(shifted, mature[1:])
                self.assertEqual(len(shifted), len(mature) - 1)

    def test_negative_left_shift_extends_the_five_prime_end(self):
        for fname in self._all_species_files():
            with self.subTest(catalog=fname):
                _, mature, pre, start, end = self._resolve(fname)
                shifted, _ = apply_shift(pre, mature, start, end, "-1|0")
                self.assertTrue(shifted.endswith(mature))
                self.assertEqual(len(shifted), len(mature) + 1)
                self.assertEqual(shifted[1:], mature)

    def test_right_shift_extends_only_the_three_prime_end(self):
        for fname in self._all_species_files():
            with self.subTest(catalog=fname):
                _, mature, pre, start, end = self._resolve(fname)
                shifted, _ = apply_shift(pre, mature, start, end, "0|1")
                self.assertTrue(shifted.startswith(mature))
                self.assertEqual(len(shifted), len(mature) + 1)

    def test_human_and_mouse_mir21_agree_despite_the_stored_coordinates(self):
        """Same miRNA, same mature sequence, opposite conventions on disk."""
        for fname, mirna_id in ((LEGACY_HUMAN_FILE, "hsa-miR-21-5p"),
                                ("mmu_mature_pre_mirna_ext.json",
                                 "mmu-miR-21a-5p")):
            with self.subTest(mirna=mirna_id):
                path = os.path.join(SHIPPED_RESOURCES_DIR, fname)
                mature, pre, start, end, _ = find_mirna_sequence(
                    load_json(path), mirna_id,
                )
                self.assertEqual(mature, MIR21)
                self.assertEqual(pre[start - 1:end], MIR21)
                self.assertEqual(apply_shift(pre, mature, start, end, "0|0")[0],
                                 MIR21)
                self.assertEqual(apply_shift(pre, mature, start, end, "1|0")[0],
                                 "AGCUUAUCAGACUGAUGUUGA")

    def test_rejects_a_shift_past_the_precursor_start(self):
        _, mature, pre, start, end = self._resolve(LEGACY_HUMAN_FILE)
        with self.assertRaises(ValueError):
            apply_shift(pre, mature, start, end, "{}|0".format(-start))

    def test_rejects_a_shift_past_the_precursor_end(self):
        _, mature, pre, start, end = self._resolve(LEGACY_HUMAN_FILE)
        with self.assertRaises(ValueError):
            apply_shift(pre, mature, start, end, "0|{}".format(len(pre)))


@unittest.skipUnless(
    os.path.isdir(SHIPPED_RESOURCES_DIR),
    "shipped resources dir not present (expected at v2/opt/resources)",
)
class CliShiftTest(unittest.TestCase):
    """End-to-end, the way app_v1/app.py invokes the script: the FASTA a 0|0
    variant produces must repeat the WT sequence."""

    SCRIPT = os.path.join(_REPO_ROOT, "v2", "mirna_processing.py")

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.env = dict(os.environ)
        self.env["ISOTAR_RESOURCES_DIR"] = SHIPPED_RESOURCES_DIR

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, mirna_id, extra=()):
        out = os.path.join(self.tmp, "mirna.fa")
        cmd = [sys.executable, self.SCRIPT, mirna_id, "-o", out] + list(extra)
        subprocess.check_call(cmd, env=self.env)
        with open(out) as f:
            return [l.strip() for l in f if l.strip()]

    def _assert_zero_shift_matches_wt(self, fname):
        mirna_id = _first_clean_id(os.path.join(SHIPPED_RESOURCES_DIR, fname))
        lines = self._run(mirna_id, ("-s", "0|0"))
        self.assertEqual(len(lines), 4,
                         "expected WT + shifted record, got {}".format(lines))
        self.assertEqual(lines[0], ">{},WT".format(mirna_id))
        self.assertEqual(lines[2], ">{},0|0,shifted".format(mirna_id))
        self.assertEqual(lines[3], lines[1],
                         "0|0 variant differs from WT for {}".format(mirna_id))

    def test_zero_shift_matches_wt_for_human(self):
        self._assert_zero_shift_matches_wt(LEGACY_HUMAN_FILE)

    def test_zero_shift_matches_wt_for_mouse(self):
        self._assert_zero_shift_matches_wt("mmu_mature_pre_mirna_ext.json")

    def test_zero_shift_matches_wt_for_worm(self):
        self._assert_zero_shift_matches_wt("cel_mature_pre_mirna_ext.json")


if __name__ == "__main__":
    unittest.main()
