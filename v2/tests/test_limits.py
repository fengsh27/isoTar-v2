"""Tests for app_v1.limits.validate_cores and env-driven defaults.

Run with:
    python3 -m unittest v2.tests.test_limits
"""

import importlib
import os
import sys
import unittest

# Make repo importable when running this file directly.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)


class ValidateCoresTests(unittest.TestCase):
    """validate_cores() takes a client-supplied value and returns (cores, error).
    error is None on success or (message, status_code) on failure."""

    def setUp(self):
        from app_v1.limits import validate_cores  # noqa: F401
        self.validate = validate_cores

    # ---- valid inputs ----

    def test_none_uses_default(self):
        cores, err = self.validate(None, max_cores=8, default_cores=8)
        self.assertEqual(cores, 8)
        self.assertIsNone(err)

    def test_int_within_range(self):
        cores, err = self.validate(4, max_cores=8, default_cores=8)
        self.assertEqual(cores, 4)
        self.assertIsNone(err)

    def test_int_string_within_range(self):
        cores, err = self.validate("4", max_cores=8, default_cores=8)
        self.assertEqual(cores, 4)
        self.assertIsNone(err)

    def test_string_with_whitespace_stripped(self):
        cores, err = self.validate("  4  ", max_cores=8, default_cores=8)
        self.assertEqual(cores, 4)
        self.assertIsNone(err)

    def test_lower_boundary(self):
        cores, err = self.validate(1, max_cores=8, default_cores=8)
        self.assertEqual(cores, 1)
        self.assertIsNone(err)

    def test_upper_boundary(self):
        cores, err = self.validate(8, max_cores=8, default_cores=8)
        self.assertEqual(cores, 8)
        self.assertIsNone(err)

    # ---- out-of-range integers ----

    def test_zero_rejected(self):
        cores, err = self.validate(0, max_cores=8, default_cores=8)
        self.assertIsNone(cores)
        self.assertEqual(err, ("cores must be between 1 and 8", 400))

    def test_negative_rejected(self):
        cores, err = self.validate(-1, max_cores=8, default_cores=8)
        self.assertIsNone(cores)
        self.assertEqual(err, ("cores must be between 1 and 8", 400))

    def test_exceeds_max_rejected(self):
        cores, err = self.validate(9, max_cores=8, default_cores=8)
        self.assertIsNone(cores)
        self.assertEqual(err, ("cores must be between 1 and 8", 400))

    def test_huge_value_rejected(self):
        cores, err = self.validate(99999, max_cores=8, default_cores=8)
        self.assertIsNone(cores)
        self.assertEqual(err[1], 400)

    # ---- type rejection ----

    def test_float_rejected(self):
        cores, err = self.validate(3.14, max_cores=8, default_cores=8)
        self.assertIsNone(cores)
        self.assertEqual(err, ("cores must be a positive integer", 400))

    def test_non_numeric_string_rejected(self):
        cores, err = self.validate("abc", max_cores=8, default_cores=8)
        self.assertIsNone(cores)
        self.assertEqual(err, ("cores must be a positive integer", 400))

    def test_empty_string_rejected(self):
        cores, err = self.validate("", max_cores=8, default_cores=8)
        self.assertIsNone(cores)
        self.assertEqual(err, ("cores must be a positive integer", 400))

    def test_bool_true_rejected(self):
        # True is technically int(1), but accepting it would be surprising.
        cores, err = self.validate(True, max_cores=8, default_cores=8)
        self.assertIsNone(cores)
        self.assertEqual(err, ("cores must be a positive integer", 400))

    def test_bool_false_rejected(self):
        cores, err = self.validate(False, max_cores=8, default_cores=8)
        self.assertIsNone(cores)
        self.assertEqual(err, ("cores must be a positive integer", 400))

    def test_list_rejected(self):
        cores, err = self.validate([4], max_cores=8, default_cores=8)
        self.assertIsNone(cores)
        self.assertEqual(err, ("cores must be a positive integer", 400))

    def test_dict_rejected(self):
        cores, err = self.validate({"cores": 4}, max_cores=8, default_cores=8)
        self.assertIsNone(cores)
        self.assertEqual(err, ("cores must be a positive integer", 400))

    # ---- max parameter respected ----

    def test_custom_max_used_in_error(self):
        cores, err = self.validate(20, max_cores=16, default_cores=8)
        self.assertIsNone(cores)
        self.assertEqual(err, ("cores must be between 1 and 16", 400))

    def test_value_at_custom_max_accepted(self):
        cores, err = self.validate(16, max_cores=16, default_cores=8)
        self.assertEqual(cores, 16)
        self.assertIsNone(err)

    # ---- module defaults applied when args omitted ----

    def test_uses_module_defaults_when_max_omitted(self):
        from app_v1 import limits as L
        cores, err = self.validate(None)
        self.assertEqual(cores, L.DEFAULT_CORES)
        self.assertIsNone(err)


class EnvDrivenDefaultsTests(unittest.TestCase):
    """Module-level constants honor env vars at import time."""

    _ENV_KEYS = (
        "ISOTAR_MAX_CORES_PER_JOB",
        "ISOTAR_DEFAULT_CORES",
        "ISOTAR_MAX_CONCURRENT_JOBS",
    )

    def setUp(self):
        # Snapshot env, will restore in tearDown.
        self._saved = {k: os.environ.get(k) for k in self._ENV_KEYS}

    def tearDown(self):
        for k, v in self._saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        # Reload to leave the cached module in a sane state for other tests.
        import app_v1.limits
        importlib.reload(app_v1.limits)

    def _reload(self):
        import app_v1.limits
        return importlib.reload(app_v1.limits)

    def test_built_in_defaults(self):
        for k in self._ENV_KEYS:
            os.environ.pop(k, None)
        L = self._reload()
        self.assertEqual(L.MAX_CORES_PER_JOB,   8)
        self.assertEqual(L.DEFAULT_CORES,       8)
        self.assertEqual(L.MAX_CONCURRENT_JOBS, 4)

    def test_env_overrides(self):
        os.environ["ISOTAR_MAX_CORES_PER_JOB"]   = "16"
        os.environ["ISOTAR_DEFAULT_CORES"]       = "12"
        os.environ["ISOTAR_MAX_CONCURRENT_JOBS"] = "2"
        L = self._reload()
        self.assertEqual(L.MAX_CORES_PER_JOB,   16)
        self.assertEqual(L.DEFAULT_CORES,       12)
        self.assertEqual(L.MAX_CONCURRENT_JOBS, 2)

    def test_invalid_env_falls_back_to_default(self):
        os.environ["ISOTAR_MAX_CORES_PER_JOB"] = "not-a-number"
        L = self._reload()
        self.assertEqual(L.MAX_CORES_PER_JOB, 8)

    def test_zero_or_negative_env_falls_back_to_default(self):
        os.environ["ISOTAR_MAX_CORES_PER_JOB"]   = "0"
        os.environ["ISOTAR_MAX_CONCURRENT_JOBS"] = "-2"
        L = self._reload()
        self.assertEqual(L.MAX_CORES_PER_JOB, 8)
        self.assertEqual(L.MAX_CONCURRENT_JOBS, 4)

    def test_empty_string_env_falls_back_to_default(self):
        os.environ["ISOTAR_MAX_CORES_PER_JOB"] = ""
        L = self._reload()
        self.assertEqual(L.MAX_CORES_PER_JOB, 8)


if __name__ == "__main__":
    unittest.main()
