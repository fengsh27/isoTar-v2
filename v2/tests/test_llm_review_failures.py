"""The model review must not report success when it read nothing.

Testers saw "read 0 papers for N queued triplets, N could not be read" under
a succeeded review, with every row saying the triplet was never queued. The
cause was python3.11 trusting no CA certificate (fixed in the Dockerfile);
the reporting bug is covered here: a hit that fails is written as a
`not_read` verdict with the reason, and a run in which nothing was read exits
non-zero so pipeline._run_llm_review marks the review failed.

llm_review.py imports the Anthropic SDK only inside main(), so its helpers are
importable under any Python 3. The pipeline test runs a stand-in script.

Needs the backend app deps for the pipeline half:
    python3 -m unittest v2.tests.test_llm_review_failures
"""

import importlib.util
import json
import os
import shutil
import sys
import tempfile
import textwrap
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

_LLM = os.path.join(_REPO, "app_v1", "verification", "llm_review.py")
_spec = importlib.util.spec_from_file_location("llm_review_under_test", _LLM)
llm_review = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(llm_review)

from app_v1.verification import pipeline  # noqa: E402


class SummarizeTest(unittest.TestCase):

    def test_clean_run(self):
        summary, rc = llm_review.summarize(3, 1, 0, 0)
        self.assertEqual(rc, 0)
        self.assertEqual(summary, {"reviewed": 3, "cached": 1, "errors": 0, "skipped_over_budget": 0})

    def test_partial_failures_still_succeed(self):
        summary, rc = llm_review.summarize(2, 0, 3, 0, first_error="fetch failed pmid=1: boom")
        self.assertEqual(rc, 0)
        self.assertNotIn("error", summary)

    def test_cache_only_run_is_not_a_failure(self):
        _, rc = llm_review.summarize(0, 4, 2, 0, first_error="x")
        self.assertEqual(rc, 0)

    def test_nothing_read_is_exit_5_with_the_first_error(self):
        summary, rc = llm_review.summarize(0, 0, 12, 0, first_error="fetch failed pmid=29484127: CERTIFICATE_VERIFY_FAILED")
        self.assertEqual(rc, 5)
        self.assertEqual(summary["error"],
                         "none of the 12 queued paper(s) could be read; first error: "
                         "fetch failed pmid=29484127: CERTIFICATE_VERIFY_FAILED")

    def test_fatal_wins(self):
        summary, rc = llm_review.summarize(0, 0, 2, 0, fatal="gateway rejected the credentials", first_error="x")
        self.assertEqual(rc, 4)
        self.assertEqual(summary["error"], "gateway rejected the credentials")

    def test_empty_queue_is_fine(self):
        _, rc = llm_review.summarize(0, 0, 0, 0)
        self.assertEqual(rc, 0)


class FailureRowTest(unittest.TestCase):

    def test_row_shape(self):
        row = llm_review.failure_row("fetch failed pmid=1: tls")
        self.assertEqual(row["verdict"], "not_read")
        self.assertEqual(row["note"], "fetch failed pmid=1: tls")
        self.assertEqual(sorted(row), sorted(c for c in llm_review.COLUMNS if c not in ("id", "scope", "pmid")))
        self.assertIn("not_read", llm_review.VERDICTS)


_FAKE_SCRIPT = textwrap.dedent('''
    import json, sys
    queue, out = sys.argv[1], sys.argv[2]
    mode = open(queue).read().strip()
    with open(out, "w") as f:
        f.write("id\\tscope\\tpmid\\tverdict\\tarm_match\\ttissue\\tevidence_type\\tsentences\\tnote\\n")
        f.write("1\\tmir_gene\\t111\\tnot_read\\tunstated\\t\\tnone\\t\\tfetch failed pmid=111: tls\\n")
        if mode == "partial":
            f.write("2\\tsame_paper\\t222\\treported_interaction\\tsame\\tlung\\treporter\\tquote\\t\\n")
    print("reviewed id=2", file=sys.stderr)
    if mode == "partial":
        print(json.dumps({"reviewed": 1, "cached": 0, "errors": 1}))
        sys.exit(0)
    print(json.dumps({"reviewed": 0, "cached": 0, "errors": 1,
                      "error": "none of the 1 queued paper(s) could be read; first error: fetch failed pmid=111: tls"}))
    sys.exit(5)
''')


class PipelineReviewStatusTest(unittest.TestCase):
    """pipeline._run_llm_review turns the script's exit code and last line into
    the review record the API serves."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        os.makedirs(pipeline.verification_dir(self.tmp))
        self.script = os.path.join(self.tmp, "fake_llm_review.py")
        with open(self.script, "w") as f:
            f.write(_FAKE_SCRIPT)
        self.logs = []

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, mode):
        with open(os.path.join(pipeline.verification_dir(self.tmp), "review_queue.json"), "w") as f:
            f.write(mode)
        config = {"python": sys.executable, "llm_script": self.script, "model": "test-model",
                  "max_papers": 5, "llm_timeout": 30}
        return pipeline._run_llm_review(self.tmp, config, "lung", lambda *a: self.logs.append(a))

    def test_nothing_read_is_a_failed_review_with_the_reason(self):
        review = self._run("none")
        self.assertEqual(review["status"], "failed")
        self.assertEqual(review["error"],
                         "none of the 1 queued paper(s) could be read; first error: fetch failed pmid=111: tls")
        self.assertEqual(review["reviewed"], 0)
        self.assertEqual(review["errors"], 1)
        # The not_read row is still served so the UI can show the reason per triplet.
        self.assertEqual([v["verdict"] for v in review["verdicts"]], ["not_read"])
        self.assertEqual(review["verdicts"][0]["note"], "fetch failed pmid=111: tls")

    def test_partial_failures_succeed_with_counts_and_rows(self):
        review = self._run("partial")
        self.assertEqual(review["status"], "succeeded")
        self.assertIsNone(review["error"])
        self.assertEqual((review["reviewed"], review["errors"]), (1, 1))
        self.assertEqual([v["verdict"] for v in review["verdicts"]], ["not_read", "reported_interaction"])

    def test_stderr_lands_in_the_log_file(self):
        self._run("none")
        with open(os.path.join(pipeline.verification_dir(self.tmp), "llm_review.log")) as f:
            self.assertIn("reviewed id=2", f.read())


if __name__ == "__main__":
    unittest.main()
