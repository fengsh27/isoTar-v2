"""Regression tests for species-aware miRNA metadata resolution.

v2/mirna_processing.py was hardcoded to the human metadata JSON
(`/opt/resources/mature_pre_mirna_ext.json`). It now resolves
`<species>_mature_pre_mirna_ext.json` from the miRBase prefix in the
miRNA id (hsa-, mmu-, dre-, ...), falling back to the un-prefixed
legacy file only for human / unparseable ids.

Runnable with:  python3 -m unittest v2.tests.test_mirna_metadata_resolver
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from v2.mirna_processing import (
    LEGACY_HUMAN_FILE,
    _species_from_mirna_id,
    find_mirna_sequence,
    load_json,
    resolve_metadata_path,
)

# Real shipped resources directory; used by the per-species integration tests.
SHIPPED_RESOURCES_DIR = os.path.join(_REPO_ROOT, "v2", "opt", "resources")

# Species we ship a metadata JSON for, beyond the un-prefixed human file.
SHIPPED_NON_HUMAN_SPECIES = (
    "cel", "cfa", "dme", "dre", "mdo", "mml", "mmu", "ptr", "rno",
)

VALID_NUC_RE = re.compile(r"^[ACGUT]+$")


class SpeciesPrefixTest(unittest.TestCase):
    def test_extracts_lowercase_three_letter_prefix(self):
        self.assertEqual(_species_from_mirna_id("hsa-miR-495-3p"), "hsa")
        self.assertEqual(_species_from_mirna_id("mmu-let-7g-5p"), "mmu")
        self.assertEqual(_species_from_mirna_id("dre-miR-10a-5p"), "dre")
        self.assertEqual(_species_from_mirna_id("HSA-miR-1"), "hsa")  # normalised

    def test_returns_none_for_unparseable_ids(self):
        self.assertIsNone(_species_from_mirna_id(""))
        self.assertIsNone(_species_from_mirna_id(None))
        self.assertIsNone(_species_from_mirna_id("no-hyphenless"))  # actually fine
        # Not 3 letters:
        self.assertIsNone(_species_from_mirna_id("ab-miR-1"))
        self.assertIsNone(_species_from_mirna_id("abcd-miR-1"))
        # Not all letters:
        self.assertIsNone(_species_from_mirna_id("hs1-miR-1"))


class ResolveMetadataPathTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        # Touch the files we want resolve_metadata_path to find.
        for name in (LEGACY_HUMAN_FILE,
                     "mmu_mature_pre_mirna_ext.json",
                     "dre_mature_pre_mirna_ext.json"):
            open(os.path.join(self.tmp, name), "w").close()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_override_wins(self):
        custom = os.path.join(self.tmp, "custom.json")
        open(custom, "w").close()
        self.assertEqual(
            resolve_metadata_path("mmu-let-7g-5p",
                                  resources_dir=self.tmp,
                                  override=custom),
            custom,
        )

    def test_species_file_used_when_present(self):
        self.assertEqual(
            resolve_metadata_path("mmu-let-7g-5p", resources_dir=self.tmp),
            os.path.join(self.tmp, "mmu_mature_pre_mirna_ext.json"),
        )
        self.assertEqual(
            resolve_metadata_path("dre-miR-10a-5p", resources_dir=self.tmp),
            os.path.join(self.tmp, "dre_mature_pre_mirna_ext.json"),
        )

    def test_human_prefers_species_file_then_legacy(self):
        # No hsa_*.json present here -> falls back to LEGACY_HUMAN_FILE.
        self.assertEqual(
            resolve_metadata_path("hsa-miR-495-3p", resources_dir=self.tmp),
            os.path.join(self.tmp, LEGACY_HUMAN_FILE),
        )
        # When hsa_*.json IS present, it wins over the legacy filename.
        hsa_file = os.path.join(self.tmp, "hsa_mature_pre_mirna_ext.json")
        open(hsa_file, "w").close()
        self.assertEqual(
            resolve_metadata_path("hsa-miR-495-3p", resources_dir=self.tmp),
            hsa_file,
        )

    def test_unknown_species_raises_not_silent_human_fallback(self):
        # 'xxx' has no file; we must NOT silently load the human file with
        # the wrong species' sequences.
        with self.assertRaises(FileNotFoundError):
            resolve_metadata_path("xxx-miR-1", resources_dir=self.tmp)

    def test_unparseable_id_falls_back_to_legacy(self):
        # Preserves prior behavior for malformed ids: use the legacy file.
        self.assertEqual(
            resolve_metadata_path("garbage", resources_dir=self.tmp),
            os.path.join(self.tmp, LEGACY_HUMAN_FILE),
        )

    def test_missing_legacy_for_human_raises(self):
        os.remove(os.path.join(self.tmp, LEGACY_HUMAN_FILE))
        with self.assertRaises(FileNotFoundError):
            resolve_metadata_path("hsa-miR-1", resources_dir=self.tmp)


def _first_id(json_path):
    """Return the first miRNA id key in a metadata file."""
    with open(json_path, "r") as f:
        data = json.load(f)
    return next(iter(data))


def _first_clean_id(json_path):
    """Return the first id in `json_path` whose entry has a single precursor
    with non-empty mature/extended-precursor sequences.

    Used by end-to-end tests so a few entries with empty `ext_pre_seq`
    (separate data-quality issue in the shipped JSONs, e.g. cfa-miR-448,
    ptr-miR-15b) and ids with multiple precursors that need --pre-id
    (e.g. rno-miR-322-5p) don't mask the resolver-level wiring we're
    actually testing here.
    """
    with open(json_path, "r") as f:
        data = json.load(f)
    for mirna_id, entries in data.items():
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


@unittest.skipUnless(
    os.path.isdir(SHIPPED_RESOURCES_DIR),
    "shipped resources dir not present (expected at v2/opt/resources)",
)
class ShippedSpeciesFilesTest(unittest.TestCase):
    """End-to-end checks against the real metadata files shipped under
    v2/opt/resources/. These guard against accidental renames, mis-prefixed
    ids inside a species file, or species being dropped from the resolver.
    """

    def test_every_shipped_species_file_resolves_correctly(self):
        for species in SHIPPED_NON_HUMAN_SPECIES:
            with self.subTest(species=species):
                expected = os.path.join(
                    SHIPPED_RESOURCES_DIR,
                    "{}_mature_pre_mirna_ext.json".format(species),
                )
                self.assertTrue(
                    os.path.exists(expected),
                    "missing shipped metadata file: {}".format(expected),
                )

                # First id in the file should start with the species prefix.
                mirna_id = _first_id(expected)
                self.assertEqual(
                    _species_from_mirna_id(mirna_id), species,
                    "first id in {} does not start with prefix '{}-': {}".format(
                        os.path.basename(expected), species, mirna_id,
                    ),
                )

                # The resolver picks this file from just the miRNA id.
                self.assertEqual(
                    resolve_metadata_path(mirna_id,
                                          resources_dir=SHIPPED_RESOURCES_DIR),
                    expected,
                )

    def test_human_id_resolves_to_legacy_file_in_shipped_dir(self):
        legacy = os.path.join(SHIPPED_RESOURCES_DIR, LEGACY_HUMAN_FILE)
        self.assertTrue(os.path.exists(legacy),
                        "missing shipped human metadata: {}".format(legacy))
        mirna_id = _first_id(legacy)
        self.assertEqual(_species_from_mirna_id(mirna_id), "hsa")
        self.assertEqual(
            resolve_metadata_path(mirna_id, resources_dir=SHIPPED_RESOURCES_DIR),
            legacy,
        )

    def test_each_species_yields_a_valid_mature_sequence(self):
        """Round-trip: resolver -> load_json -> find_mirna_sequence.

        Picks a clean (single-precursor, populated) id per species so this
        test exercises every shipped species end-to-end without tripping on
        unrelated data-quality issues in specific entries.
        """
        species_files = [(LEGACY_HUMAN_FILE, "hsa")] + [
            ("{}_mature_pre_mirna_ext.json".format(s), s)
            for s in SHIPPED_NON_HUMAN_SPECIES
        ]
        for fname, species in species_files:
            with self.subTest(species=species):
                json_path = os.path.join(SHIPPED_RESOURCES_DIR, fname)
                mirna_id = _first_clean_id(json_path)
                # Sanity: the id we picked still maps to the right species.
                self.assertEqual(_species_from_mirna_id(mirna_id), species)

                resolved = resolve_metadata_path(
                    mirna_id, resources_dir=SHIPPED_RESOURCES_DIR,
                )
                data = load_json(resolved)
                mature, pre, start, end, _ = find_mirna_sequence(data, mirna_id)

                self.assertIsNotNone(mature,
                                     "no mature sequence for {}".format(mirna_id))
                self.assertTrue(VALID_NUC_RE.match(mature),
                                "invalid nucleotides in mature seq for {}: {!r}"
                                .format(mirna_id, mature))
                self.assertTrue(VALID_NUC_RE.match(pre),
                                "invalid nucleotides in precursor for {}".format(mirna_id))
                self.assertGreaterEqual(start, 0)
                self.assertGreater(end, start)
                self.assertLessEqual(end, len(pre))


@unittest.skipUnless(
    os.path.isdir(SHIPPED_RESOURCES_DIR),
    "shipped resources dir not present (expected at v2/opt/resources)",
)
class MirnaProcessingCliPerSpeciesTest(unittest.TestCase):
    """Invoke v2/mirna_processing.py end-to-end (the same way app_v1/app.py
    calls it) for a representative non-human species and verify the FASTA
    output. This is the closest thing to a real job submission without
    standing up Celery / the API."""

    SCRIPT = os.path.join(_REPO_ROOT, "v2", "mirna_processing.py")

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.env = dict(os.environ)
        self.env["ISOTAR_RESOURCES_DIR"] = SHIPPED_RESOURCES_DIR

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _clean_id(self, fname):
        return _first_clean_id(os.path.join(SHIPPED_RESOURCES_DIR, fname))

    def _run(self, mirna_id, extra=()):
        out = os.path.join(self.tmp, "mirna.fa")
        cmd = [sys.executable, self.SCRIPT, mirna_id, "-o", out] + list(extra)
        subprocess.check_call(cmd, env=self.env)
        with open(out) as f:
            return f.read()

    def _assert_fasta_wt(self, content, mirna_id):
        lines = [l for l in content.splitlines() if l.strip()]
        self.assertEqual(len(lines), 2,
                         "expected 1 header + 1 seq, got: {}".format(content))
        self.assertEqual(lines[0], ">{},WT".format(mirna_id))
        self.assertTrue(VALID_NUC_RE.match(lines[1]),
                        "invalid nucleotides in FASTA seq: {!r}".format(lines[1]))

    def test_runs_for_mouse(self):
        mirna_id = self._clean_id("mmu_mature_pre_mirna_ext.json")
        self._assert_fasta_wt(self._run(mirna_id), mirna_id)

    def test_runs_for_zebrafish(self):
        mirna_id = self._clean_id("dre_mature_pre_mirna_ext.json")
        self._assert_fasta_wt(self._run(mirna_id), mirna_id)

    def test_runs_for_rat(self):
        mirna_id = self._clean_id("rno_mature_pre_mirna_ext.json")
        self._assert_fasta_wt(self._run(mirna_id), mirna_id)

    def test_runs_for_worm(self):
        mirna_id = self._clean_id("cel_mature_pre_mirna_ext.json")
        self._assert_fasta_wt(self._run(mirna_id), mirna_id)

    def test_runs_for_fly(self):
        mirna_id = self._clean_id("dme_mature_pre_mirna_ext.json")
        self._assert_fasta_wt(self._run(mirna_id), mirna_id)

    def test_runs_for_human_via_legacy_file(self):
        mirna_id = _first_clean_id(
            os.path.join(SHIPPED_RESOURCES_DIR, LEGACY_HUMAN_FILE)
        )
        self._assert_fasta_wt(self._run(mirna_id), mirna_id)

    def test_metadata_override_flag_wins(self):
        """The --metadata flag should override species-based resolution."""
        # Resolve a mouse id but force the metadata path to the rat file --
        # the run should fail (mouse id not in rat data), proving --metadata
        # actually drove the lookup rather than the species prefix.
        mouse_id = self._clean_id("mmu_mature_pre_mirna_ext.json")
        rat_file = os.path.join(SHIPPED_RESOURCES_DIR,
                                "rno_mature_pre_mirna_ext.json")
        out = os.path.join(self.tmp, "mirna.fa")
        proc = subprocess.run(
            [sys.executable, self.SCRIPT, mouse_id, "-o", out,
             "--metadata", rat_file],
            env=self.env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        self.assertNotEqual(proc.returncode, 0,
                            "mouse id against rat metadata should fail")

    def test_env_var_overrides_resources_dir(self):
        """Setting ISOTAR_RESOURCES_DIR to a directory without the species
        file makes resolution fail, even though the shipped dir would work.
        Proves the env var is the only knob the caller needs."""
        empty = tempfile.mkdtemp()
        try:
            env = dict(os.environ)
            env["ISOTAR_RESOURCES_DIR"] = empty
            out = os.path.join(self.tmp, "mirna.fa")
            proc = subprocess.run(
                [sys.executable, self.SCRIPT, "mmu-let-7g-5p", "-o", out],
                env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn(b"mmu", proc.stderr + proc.stdout)
        finally:
            shutil.rmtree(empty, ignore_errors=True)

    def test_unknown_species_fails_loudly(self):
        # No xxx_*.json shipped -> script must exit non-zero, not silently
        # use the human file.
        out = os.path.join(self.tmp, "mirna.fa")
        proc = subprocess.run(
            [sys.executable, self.SCRIPT, "xxx-miR-1", "-o", out],
            env=self.env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse(os.path.exists(out),
                         "FASTA should not be written when resolver fails")
        self.assertIn(b"xxx", proc.stderr + proc.stdout)


if __name__ == "__main__":
    unittest.main()
