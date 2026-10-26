#!/usr/bin/env python3.11
"""Read the papers behind the PubMed hits and record one verdict per (triplet, scope, PMID).

Runs under python3.11 (the Anthropic SDK needs >= 3.10; the rest of app_v1 is 3.6).
Invoked by app_v1/verification/pipeline.py; usable by hand:

    python3.11 llm_review.py review_queue.json verdicts.tsv --model <name> [--max-papers 40]
                            [--cache verdicts.db] [--tissue lung]

Environment:
    ANTHROPIC_BASE_URL / ANTHROPIC_AUTH_TOKEN   the LiteLLM (or other) gateway and its key,
                                               read by the SDK itself
    NCBI_API_KEY / NCBI_EMAIL                  for efetch (abstract and PMC full text)

For each queued hit the abstract plus open-access full text is fetched, the
sentences that contain the symbols are extracted, and the model classifies the
hit with the rubric below, returning a structured verdict that always carries
the quoted sentence deciding it. Verdicts are cached by (PMID, scope, symbols,
model, rubric version) so a paper read for one job is free for the next.

Exit codes: 0 ok, 2 bad arguments, 3 SDK not installed, 4 gateway/auth failure.
The last stdout line is a JSON summary: {"reviewed": n, "cached": c, "errors": e}.
"""
import argparse
import json
import os
import sqlite3
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import papers  # noqa: E402  (stdlib sibling module)

RUBRIC_VERSION = 1
VERDICTS = ("reported_interaction", "opposite_arm", "co_mention", "look_alike", "editorial_notice", "not_accessible")
ARM_MATCH = ("same", "opposite", "unstated", "na")
EVIDENCE = ("reporter", "rip", "expression", "computational", "review", "none")
COLUMNS = ["id", "scope", "pmid", "verdict", "arm_match", "tissue", "evidence_type", "sentences", "note"]

SYSTEM_PROMPT = """You classify PubMed hits for predicted ceRNA triplets (miRNA -> mRNA -> lncRNA).
A keyword search found that two or three symbols co-occur in a paper. Co-occurrence is not a reported
interaction. Decide, from the sentences you are given, what the co-occurrence is.

Scope of a hit:
- same_paper: the whole axis (miRNA family, gene and lncRNA) must be reported together.
- mir_gene: the miRNA and the gene. mir_lnc: the miRNA and the lncRNA. gene_lnc: the gene and the lncRNA.

Procedure:
1. Find the sentence(s) in which the two (or three) members of the scope are linked by a verb: targets,
   sponges, binds, regulates, correlates, is a target of, axis. If no sentence links them, the paper is a
   co-mention (gene lists, signatures, network figures, two separate findings).
2. Check the mature arm. The search is family-level, so a paper on miR-141-3p satisfies a miR-141-5p
   prediction. If the paper names the opposite arm, the verdict is opposite_arm even when the mechanism is
   otherwise the predicted one. If the arm is never stated, arm_match is "unstated". For gene_lnc, "na".
3. Check the direction when stated (a miRNA represses its target; a sponging lncRNA raises the target) and
   record the tissue as the paper names it (for example "gastric cancer", "NSCLC", "cell lines only").
4. Record the evidence type from the sentence: reporter (luciferase), rip (RNA immunoprecipitation,
   pull-down, MS2), expression (correlation or knockdown effect only), computational (predicted or
   inferred bioinformatically), review (secondary citation), none.

Verdict classes:
- reported_interaction: the paper reports the predicted relationship between the members of the scope
  (for same_paper, the complete axis with this miRNA family, gene and lncRNA).
- opposite_arm: the relationship is reported for the other mature arm.
- co_mention: both symbols appear but nothing links them.
- look_alike: the symbol matched inside another token (ZEB2 inside ZEB2-AS1, KIT in "Cell Counting Kit-8",
  MET in "metastasis"), or a natural antisense pair that always co-occurs without a reported regulation.
- editorial_notice: a retraction, expression of concern or erratum that lists titles.
- not_accessible: the abstract alone cannot decide and the full text was not available; say what the
  abstract does state.

A lung ceRNA axis of the same lncRNA and gene through a different miRNA is NOT a same_paper match; it is a
gene_lnc reported_interaction with a note naming the published miRNA.

Always quote verbatim, in "sentences", the sentence or sentences that decide the verdict, joined by " | ",
trimmed to what is needed. A verdict without a quote is worth nothing: a reader confirms your verdict from
the quote without re-reading the paper. Put anything else the reader needs in "note"."""


def _try_import_sdk():
    try:
        import anthropic  # noqa: F401
        from pydantic import BaseModel  # noqa: F401
    except ImportError as e:
        print(json.dumps({"reviewed": 0, "cached": 0, "errors": 0,
                          "error": "python3.11 is missing the anthropic/pydantic packages: %s" % e}))
        sys.exit(3)


class VerdictCache(object):
    def __init__(self, path):
        self.conn = sqlite3.connect(path, timeout=60) if path else sqlite3.connect(":memory:")
        self.conn.execute("CREATE TABLE IF NOT EXISTS verdicts (key TEXT PRIMARY KEY, value TEXT NOT NULL, created_at INTEGER NOT NULL)")
        self.conn.commit()

    @staticmethod
    def key(pmid, scope, mirna_family, gene, lncrna, model):
        return "|".join([str(pmid), scope, mirna_family, gene, lncrna or "", model, "r%d" % RUBRIC_VERSION])

    def get(self, key):
        row = self.conn.execute("SELECT value FROM verdicts WHERE key = ?", (key,)).fetchone()
        return json.loads(row[0]) if row else None

    def put(self, key, value):
        self.conn.execute("INSERT OR REPLACE INTO verdicts (key, value, created_at) VALUES (?, ?, ?)",
                          (key, json.dumps(value), int(time.time())))
        self.conn.commit()


def _scope_members(scope, item):
    fam, gene, lnc = item["terms"][0], item["gene"], item["lncRNA"]
    if scope == "same_paper":
        return [fam, gene, lnc]
    if scope == "mir_gene":
        return [fam, gene]
    if scope == "mir_lnc":
        return [fam, lnc]
    return [gene, lnc]


def _user_message(item, scope, paper, found, tissue_label):
    members = _scope_members(scope, item)
    lines = [
        "Predicted triplet (rank %s): %s -> %s -> %s" % (item["id"], item["miRNA"], item["gene"], item["lncRNA"]),
        "Scope to judge: %s (%s)" % (scope, " and ".join(members)),
        "Cohort tissue context: %s" % tissue_label,
        "",
        "PMID %s | %s %s" % (paper["pmid"], paper["journal"], paper["year"]),
        "Title: %s" % paper["title"],
        "",
        "ABSTRACT",
        paper["abstract"] or "(no abstract)",
        "",
        "FULL TEXT: %s" % (("%d words from %s" % (len(paper["body"].split()), paper["source"])) if paper["body"] else "not available"),
    ]
    for term, sents in found.items():
        lines.append("")
        lines.append("Sentences containing %s (%d):" % (term, len(sents)))
        for s in sents:
            lines.append("  - " + s)
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("queue")
    ap.add_argument("out")
    ap.add_argument("--model", required=True)
    ap.add_argument("--max-papers", type=int, default=40)
    ap.add_argument("--cache", default=None)
    ap.add_argument("--tissue", default="lung")
    ap.add_argument("--max-tokens", type=int, default=2048)
    a = ap.parse_args()
    if not a.model:
        print(json.dumps({"reviewed": 0, "cached": 0, "errors": 0, "error": "no model configured"}))
        sys.exit(2)

    _try_import_sdk()
    import anthropic
    from pydantic import BaseModel
    from typing import Literal

    class Verdict(BaseModel):
        verdict: Literal["reported_interaction", "opposite_arm", "co_mention", "look_alike", "editorial_notice", "not_accessible"]
        arm_match: Literal["same", "opposite", "unstated", "na"]
        tissue: str
        evidence_type: Literal["reporter", "rip", "expression", "computational", "review", "none"]
        sentences: str
        note: str

    with open(a.queue) as f:
        queue = json.load(f)
    api_key, email = os.environ.get("NCBI_API_KEY", ""), os.environ.get("NCBI_EMAIL", "")
    cache = VerdictCache(a.cache)
    client = anthropic.Anthropic()  # base URL, key/token and retries from the environment

    import csv
    exists = os.path.exists(a.out)
    outf = open(a.out, "a", newline="", encoding="utf-8")
    writer = csv.DictWriter(outf, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    if not exists:
        writer.writeheader()

    reviewed = cached = errors = skipped = 0
    papers_seen = {}
    budget = a.max_papers
    fatal = None
    for item in queue:
        for scope, pmid in item["pmids"]:
            key = VerdictCache.key(pmid, scope, item["terms"][0], item["gene"], item["lncRNA"], a.model)
            hit = cache.get(key)
            if hit is not None:
                writer.writerow(dict(hit, id=item["id"], scope=scope, pmid=pmid))
                outf.flush()
                cached += 1
                continue
            if budget <= 0:
                skipped += 1
                continue
            try:
                if pmid not in papers_seen:
                    papers_seen[pmid] = papers.fetch_paper(pmid, api_key, email)
                paper = papers_seen[pmid]
            except Exception as e:  # noqa: BLE001 -- NCBI hiccup on one paper
                errors += 1
                print("fetch failed pmid=%s: %s" % (pmid, str(e)[:200]), file=sys.stderr)
                continue
            if not paper.get("found"):
                row = {"verdict": "not_accessible", "arm_match": "unstated", "tissue": "", "evidence_type": "none",
                       "sentences": "", "note": "PMID not found in PubMed"}
                writer.writerow(dict(row, id=item["id"], scope=scope, pmid=pmid))
                outf.flush()
                continue
            found = papers.sentences_with(paper, _scope_members(scope, item))
            budget -= 1
            try:
                response = client.messages.parse(
                    model=a.model,
                    max_tokens=a.max_tokens,
                    system=SYSTEM_PROMPT,
                    messages=[{"role": "user", "content": _user_message(item, scope, paper, found, a.tissue)}],
                    output_format=Verdict,
                )
            except anthropic.AuthenticationError as e:
                fatal = "gateway rejected the credentials: %s" % e
                break
            except anthropic.NotFoundError as e:
                fatal = "model %s not found at the gateway: %s" % (a.model, e)
                break
            except (anthropic.APIStatusError, anthropic.APIConnectionError) as e:
                errors += 1
                print("model call failed pmid=%s scope=%s: %s" % (pmid, scope, str(e)[:300]), file=sys.stderr)
                continue
            if response.stop_reason == "refusal" or response.parsed_output is None:
                row = {"verdict": "not_accessible", "arm_match": "unstated", "tissue": "", "evidence_type": "none",
                       "sentences": "", "note": "model returned no verdict (stop_reason=%s)" % response.stop_reason}
            else:
                v = response.parsed_output
                row = {"verdict": v.verdict, "arm_match": v.arm_match, "tissue": v.tissue,
                       "evidence_type": v.evidence_type, "sentences": v.sentences, "note": v.note}
            cache.put(key, row)
            writer.writerow(dict(row, id=item["id"], scope=scope, pmid=pmid))
            outf.flush()
            reviewed += 1
            print("reviewed id=%s scope=%s pmid=%s verdict=%s" % (item["id"], scope, pmid, row["verdict"]), file=sys.stderr)
        if fatal:
            break
    outf.close()
    summary = {"reviewed": reviewed, "cached": cached, "errors": errors, "skipped_over_budget": skipped}
    if fatal:
        summary["error"] = fatal
        print(json.dumps(summary))
        sys.exit(4)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
