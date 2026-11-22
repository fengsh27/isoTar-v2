"""llm_review.ask_verdict: a verdict through a validated tool call.

The institutional gateway (LiteLLM -> Vertex AI) rejects structured outputs
under an organization policy, and a forced tool_choice is not supported for the
model, so the verdict tool is offered with tool_choice "auto" and its input is
validated against the Verdict model, with the exact problem sent back.

Needs pydantic (python3.11 in the image has it); skipped otherwise.

Run with:
    python3 -m unittest v2.tests.test_llm_review_ask_verdict
"""

import importlib.util
import os
import sys
import unittest
from types import SimpleNamespace

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
_LLM = os.path.join(_REPO, "app_v1", "verification", "llm_review.py")

try:
    from typing import Literal

    from pydantic import BaseModel
except ImportError:  # pragma: no cover - pydantic only lives in the python3.11 image
    BaseModel = None

if BaseModel is not None:
    _spec = importlib.util.spec_from_file_location("llm_review_ask_verdict", _LLM)
    llm_review = importlib.util.module_from_spec(_spec)
    sys.path.insert(0, os.path.dirname(_LLM))  # llm_review imports its sibling `papers`
    _spec.loader.exec_module(llm_review)

    class Verdict(BaseModel):
        verdict: Literal["reported_interaction", "co_mention", "not_accessible"]
        sentences: str
        note: str

GOOD = {"verdict": "reported_interaction", "sentences": "miR-1 targets TP53.", "note": ""}


def tool_use(payload, id_="tu_1"):
    return SimpleNamespace(type="tool_use", name="record_verdict", input=payload, id=id_)


def text(t):
    return SimpleNamespace(type="text", text=t)


class FakeClient(object):
    """Replays scripted responses and records every request."""

    def __init__(self, *responses):
        self._responses = list(responses)
        self.calls = []
        self.messages = self

    def create(self, **kwargs):
        self.calls.append(kwargs)
        stop, content = self._responses.pop(0)
        return SimpleNamespace(stop_reason=stop, content=content)


@unittest.skipIf(BaseModel is None, "pydantic not installed (it is in the python3.11 image)")
class AskVerdictTest(unittest.TestCase):
    def ask(self, client):
        return llm_review.ask_verdict(client, "claude-opus-5-5", 2048, "SYSTEM", "USER", Verdict)

    def test_valid_first_answer(self):
        client = FakeClient(("tool_use", [tool_use(GOOD)]))
        verdict, problem = self.ask(client)
        self.assertIsNone(problem)
        self.assertEqual(verdict.verdict, "reported_interaction")
        call = client.calls[0]
        # Never structured outputs or a forced choice: both are rejected by the gateway.
        self.assertNotIn("output_format", call)
        self.assertEqual(call["tool_choice"], {"type": "auto"})
        self.assertEqual(call["tools"][0]["name"], "record_verdict")
        self.assertIn("verdict", call["tools"][0]["input_schema"]["properties"])
        self.assertTrue(call["system"].startswith("SYSTEM"))

    def test_missing_field_is_sent_back_and_fixed(self):
        client = FakeClient(("tool_use", [tool_use({"verdict": "co_mention"}, "tu_bad")]),
                            ("tool_use", [tool_use(GOOD, "tu_good")]))
        verdict, problem = self.ask(client)
        self.assertIsNone(problem)
        self.assertEqual(verdict.verdict, "reported_interaction")
        retry = client.calls[1]["messages"]
        result = retry[-1]["content"][0]
        self.assertEqual(result["type"], "tool_result")
        self.assertEqual(result["tool_use_id"], "tu_bad")
        self.assertTrue(result["is_error"])
        self.assertIn("sentences", result["content"])  # names the missing field

    def test_disallowed_value_is_rejected(self):
        client = FakeClient(("tool_use", [tool_use(dict(GOOD, verdict="direct_target"))]),
                            ("tool_use", [tool_use(GOOD)]))
        verdict, problem = self.ask(client)
        self.assertIsNone(problem)
        self.assertIn("verdict", client.calls[1]["messages"][-1]["content"][0]["content"])

    def test_plain_text_answer_is_asked_again(self):
        client = FakeClient(("end_turn", [text("It is a reported interaction.")]),
                            ("tool_use", [tool_use(GOOD)]))
        verdict, problem = self.ask(client)
        self.assertIsNone(problem)
        self.assertIn("record_verdict", client.calls[1]["messages"][-1]["content"])

    def test_gives_up_with_the_last_problem(self):
        bad = ("tool_use", [tool_use({"verdict": "co_mention"})])
        client = FakeClient(bad, bad, bad)
        verdict, problem = self.ask(client)
        self.assertIsNone(verdict)
        self.assertEqual(len(client.calls), 3)
        self.assertIn("no valid verdict after 3 attempt(s)", problem)
        self.assertIn("sentences", problem)

    def test_refusal_stops_at_once(self):
        client = FakeClient(("refusal", []))
        verdict, problem = self.ask(client)
        self.assertIsNone(verdict)
        self.assertEqual(len(client.calls), 1)
        self.assertIn("refused", problem)


if __name__ == "__main__":
    unittest.main()
