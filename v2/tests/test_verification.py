"""Tests for the literature verification package (app_v1.verification).

Pure parts only: query terms, the client's caching/retry against a fake
fetcher, roll-up, review queue and triplet building. No network, no Flask.

Run with:
    python3 -m unittest v2.tests.test_verification
"""

import json
import os
import sys
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app_v1.verification import pubmed  # noqa: E402
from app_v1.verification.pipeline import build_triplets  # noqa: E402
from app_v1.verification.review import build_review_queue  # noqa: E402
from app_v1.verification.rollup import COLUMNS, best_tier, rollup_rows, summarize  # noqa: E402
from app_v1.verification.symbols import LncrnaSymbols  # noqa: E402


class TermTests(unittest.TestCase):
    def test_family_strips_species_arm_and_variant(self):
        self.assertEqual(pubmed.mirna_family("hsa-miR-141-5p"), "miR-141")
        self.assertEqual(pubmed.mirna_family("hsa-let-7c-5p"), "let-7c")
        self.assertEqual(pubmed.mirna_family("hsa-miR-21-5p,8:A>U,modified"), "miR-21")

    def test_mirna_term_forms(self):
        t = pubmed.mirna_term("hsa-miR-141-5p")
        for form in ('"miR-141"[tiab]', '"microRNA-141"[tiab]', '"hsa-miR-141"[tiab]', '"miR 141"[tiab]'):
            self.assertIn(form, t)
        self.assertNotIn("5p", t)
        self.assertEqual(pubmed.mirna_term("hsa-let-7c-5p"), '("let-7c"[tiab] OR "hsa-let-7c"[tiab])')

    def test_tier_queries_nest_and_exclude_notices(self):
        qs = dict(pubmed.tier_queries('"A"[tiab]', '"B"[tiab]', pubmed.TISSUE["lung"]))
        self.assertIn('"lung cancer"[tiab]', qs["1"])
        self.assertIn("cancer*[tiab]", qs["2"])
        self.assertTrue(qs["3"].startswith('"A"[tiab] AND "B"[tiab]'))
        for q in qs.values():
            self.assertIn('NOT ("expression of concern"[pt]', q)

    def test_tier1_not_run_without_tissue(self):
        qs = dict(pubmed.tier_queries("x", "y", None))
        self.assertIsNone(qs["1"])
        self.assertIsNotNone(qs["2"])

    def test_resolve_tissue(self):
        self.assertEqual(pubmed.resolve_tissue("Breast")[1], "breast")
        self.assertEqual(pubmed.resolve_tissue("none"), (None, "none"))
        self.assertEqual(pubmed.resolve_tissue('"glioma"[tiab]')[1], "custom")
        with self.assertRaises(ValueError):
            pubmed.resolve_tissue("kidney")

    def test_unsearchable(self):
        self.assertTrue(pubmed.is_unsearchable("ENSG00000234456"))
        self.assertTrue(pubmed.is_unsearchable("ENST00000616796"))
        self.assertTrue(pubmed.is_unsearchable(""))
        self.assertFalse(pubmed.is_unsearchable("MAGI2-AS3"))


class FakeFetch(object):
    """Answers esearch with a count keyed by substrings of the term, and
    esummary with a fixed title. Counts calls so caching can be asserted."""

    def __init__(self, counts, fail_first=0):
        self.counts = counts
        self.calls = 0
        self.fail_first = fail_first

    def __call__(self, url):
        self.calls += 1
        if self.fail_first > 0:
            self.fail_first -= 1
            raise OSError("boom")
        if "esearch.fcgi" in url:
            count = 0
            for needle, n in self.counts.items():
                if needle in url:
                    count = n
                    break
            ids = ["%d" % (1000 + i) for i in range(min(count, 3))]
            return json.dumps({"esearchresult": {"count": str(count), "idlist": ids}})
        if "esummary.fcgi" in url:
            from urllib.parse import parse_qs, urlparse
            ids = parse_qs(urlparse(url).query)["id"][0].split(",")
            return json.dumps({"result": dict((i, {"title": "T" + i, "pubdate": "2020 Jan", "source": "J"}) for i in ids)})
        raise AssertionError("unexpected url")


def _client(counts, **kw):
    fetch = FakeFetch(counts, **kw)
    cl = pubmed.PubMedClient("k", "e", cache=pubmed.PubMedCache(), fetch=fetch, min_interval=0, max_sleep=0)
    return cl, fetch


class ClientTests(unittest.TestCase):
    def test_cache_hit_avoids_second_fetch(self):
        cl, fetch = _client({"ZEB1": 5})
        self.assertEqual(cl.esearch('"ZEB1"[tiab]')[0], 5)
        self.assertEqual(cl.esearch('"ZEB1"[tiab]')[0], 5)
        self.assertEqual(fetch.calls, 1)
        self.assertEqual(cl.live_calls, 1)

    def test_retries_then_succeeds(self):
        cl, fetch = _client({"ZEB1": 2}, fail_first=2)
        count, ids = cl.esearch('"ZEB1"[tiab]')
        self.assertEqual((count, ids), (2, ["1000", "1001"]))
        self.assertEqual(fetch.calls, 3)

    def test_gives_up_without_leaking_key(self):
        cl, fetch = _client({}, fail_first=99)
        cl.retries = 2
        with self.assertRaises(pubmed.PubMedError) as ctx:
            cl.esearch('"ZEB1"[tiab]')
        self.assertNotIn("api_key", str(ctx.exception))
        self.assertNotIn("k&", str(ctx.exception))

    def test_sqlite_cache_persists(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "c.db")
            c1 = pubmed.PubMedCache(path)
            c1.put("q", {"a": 1})
            c1.close()
            c2 = pubmed.PubMedCache(path)
            self.assertEqual(c2.get("q"), {"a": 1})
            self.assertIsNone(c2.get("missing"))
            c2.close()


TRIPLETS = [
    {"id": 1, "miRNA": "hsa-miR-200b-3p", "gene": "ZEB1", "lncRNA": "MAGI2-AS3"},
    {"id": 2, "miRNA": "hsa-miR-30a-5p", "gene": "FBXO45", "lncRNA": "ENST00000000001"},
    {"id": 3, "miRNA": "hsa-miR-205-5p", "gene": "AXIN2", "lncRNA": "PVT1"},
]


class SearchAndRollupTests(unittest.TestCase):
    def setUp(self):
        # Counts by the first matching substring of the encoded term, in order:
        # every query naming MAGI2-AS3 (its two arms and the same-paper query)
        # answers 1 record; the miR-200b/ZEB1 arm answers 40; FBXO45 3; AXIN2 0;
        # anything else (the miR-205/PVT1 arm) 0.
        self.cl, self.fetch = _client({
            "MAGI2-AS3": 1,
            "ZEB1": 40,
            "FBXO45": 3,
            "AXIN2": 0,
        })

    def test_distinct_arms_skip_unsearchable_lncrna(self):
        arms = pubmed.distinct_arms(TRIPLETS)
        kinds = [a[0] for a in arms]
        self.assertEqual(kinds.count("mir_gene"), 3)
        self.assertEqual(kinds.count("mir_lnc"), 2)
        self.assertEqual(kinds.count("gene_lnc"), 2)

    def test_rollup_rows(self):
        arms = pubmed.search_arms(self.cl, TRIPLETS, pubmed.TISSUE["lung"])
        same = pubmed.search_same_paper(self.cl, TRIPLETS)
        titles = self.cl.esummary(["1000"])
        rows = rollup_rows(TRIPLETS, arms, same, titles)
        self.assertEqual([r["id"] for r in rows], [1, 2, 3])
        r1, r2, r3 = rows
        # Triplet 1: every arm hits at tier 1 -> tier1, three arms, same-paper hit.
        self.assertEqual(r1["triplet_tier"], "tier1")
        self.assertEqual(r1["arms_supported"], 3)
        self.assertEqual(r1["same_paper_hit"], 1)
        self.assertEqual(r1["same_paper_url"], "https://pubmed.ncbi.nlm.nih.gov/1000/")
        self.assertEqual(r1["mir_gene_note"], "n=40; 2020 J: T1000")
        self.assertEqual(r1["mir_gene_pair"], "miR-200b-3p to ZEB1")
        # Triplet 2: lncRNA has no symbol -> two arms unsearchable, same-paper blank.
        self.assertEqual(r2["mir_lnc_tier"], "unsearchable")
        self.assertEqual(r2["gene_lnc_tier"], "unsearchable")
        self.assertEqual(r2["same_paper_hit"], "")
        self.assertEqual(r2["triplet_tier"], "partial")
        self.assertEqual(r2["arms_supported"], 1)
        # Triplet 3: nothing for AXIN2 arms, PVT1 arm untouched by counts -> none.
        self.assertEqual(r3["mir_gene_tier"], "none")
        self.assertEqual(r3["triplet_tier"], "none")
        self.assertEqual(r3["same_paper_hit"], 0)
        for r in rows:
            for c in COLUMNS:
                self.assertIn(c, r)
        s = summarize(rows)
        self.assertEqual((s["with_any_arm"], s["with_three_arms"], s["same_paper_hits"], s["unsearchable_lncrna"]),
                         (2, 1, 1, 1))

    def test_best_tier_prefers_lowest_number(self):
        self.assertEqual(best_tier({"1": (0, []), "2": (4, ["9"]), "3": (9, ["9", "8"])})[0], "2")
        self.assertEqual(best_tier({"1": ("not_run", []), "2": (0, []), "3": (0, [])})[0], "not_run")
        self.assertEqual(best_tier({"1": (0, []), "2": (0, []), "3": (0, [])})[0], "none")


class ReviewQueueTests(unittest.TestCase):
    def _row(self, **kw):
        base = {"id": 1, "miRNA": "hsa-miR-1-3p", "gene": "G", "lncRNA": "L", "triplet_tier": "partial",
                "arms_supported": 1, "same_paper_hit": 0, "same_paper_url": ""}
        for k in ("mir_gene", "mir_lnc", "gene_lnc"):
            base[k + "_tier"], base[k + "_url"], base[k + "_note"] = "none", "", ""
        base.update(kw)
        return base

    def test_selects_same_paper_three_arm_and_experimental_two_arm(self):
        rows = [
            self._row(id=1, same_paper_hit=1, same_paper_url="https://pubmed.ncbi.nlm.nih.gov/11/",
                      mir_gene_tier="3", mir_gene_url="https://pubmed.ncbi.nlm.nih.gov/12/"),
            self._row(id=2, arms_supported=3, triplet_tier="tier2", mir_gene_tier="2", mir_lnc_tier="2",
                      gene_lnc_tier="2", mir_gene_url="https://pubmed.ncbi.nlm.nih.gov/21/"),
            self._row(id=3, arms_supported=2, mir_gene_tier="3", mir_lnc_tier="3",
                      mir_gene_note="n=1; 2019 J: luciferase reporter assay"),
            self._row(id=4, arms_supported=2, mir_gene_tier="3", mir_lnc_tier="3", mir_gene_note="n=1; 2019 J: signature"),
            self._row(id=5, arms_supported=1, mir_gene_tier="1"),
        ]
        q = build_review_queue(rows)
        self.assertEqual([x["id"] for x in q], [1, 2, 3])
        self.assertEqual(q[0]["reason"], "same_paper")
        self.assertEqual(q[0]["pmids"], [["same_paper", "11"], ["mir_gene", "12"]])
        self.assertEqual(q[2]["reason"], "two_arms_experimental")
        self.assertEqual(q[0]["terms"], ["miR-1", "G", "L"])

    def test_two_arm_tissue_tier(self):
        rows = [self._row(id=9, arms_supported=2, mir_gene_tier="1", gene_lnc_tier="3")]
        self.assertEqual(build_review_queue(rows)[0]["reason"], "two_arms_tissue")
        # With tier 1 not run (tissue "none") a plain two-arm row is not queued.
        self.assertEqual(build_review_queue(rows, tissue_tier=""), [])


class StubSymbols(object):
    def __init__(self, table):
        self.table = table

    def symbol_for(self, lncrna_id):
        return self.table.get(lncrna_id)


class BuildTripletsTests(unittest.TestCase):
    def test_one_triplet_per_pair_and_base_mirna(self):
        network = {"pairs": [
            {"gene": "NM_1", "gene_label": "ZEB1", "lncrna": "ENST1", "score": 0.5,
             "bridges": [{"mirna": "hsa-miR-200b-3p", "score": 0.5, "is_variant": False},
                         {"mirna": "hsa-miR-200b-3p,8:A>U,modified", "score": 0.7, "is_variant": True},
                         {"mirna": "hsa-miR-141-3p", "score": 0.2, "is_variant": False}]},
            {"gene": "NM_2", "gene_label": "AXIN2", "lncrna": "ENST2", "score": None,
             "bridges": [{"mirna": "hsa-miR-205-5p,-7|1,shifted", "score": 0.1, "is_variant": True}]},
        ]}
        trip = build_triplets(network, StubSymbols({"ENST1": "MAGI2-AS3"}))
        self.assertEqual([(t["id"], t["miRNA"], t["gene"], t["lncRNA"]) for t in trip], [
            (1, "hsa-miR-200b-3p", "ZEB1", "MAGI2-AS3"),
            (2, "hsa-miR-141-3p", "ZEB1", "MAGI2-AS3"),
            (3, "hsa-miR-205-5p", "AXIN2", "ENST2"),
        ])
        # The variant outscored WT, so it is the recorded bridge for the collapsed triplet.
        self.assertEqual(trip[0]["mirna_node"], "hsa-miR-200b-3p,8:A>U,modified")
        self.assertEqual(trip[0]["bridge_is_variant"], 1)
        self.assertEqual(trip[0]["bridge_score"], 0.7)
        self.assertEqual(trip[0]["pair_score"], 0.5)
        self.assertEqual(trip[1]["bridge_is_variant"], 0)
        self.assertEqual(trip[2]["pair_score"], "")
        self.assertEqual(trip[2]["lncrna_id"], "ENST2")


class SymbolsTests(unittest.TestCase):
    def test_symbol_lookup_from_tsv_and_reference_db(self):
        import sqlite3
        with tempfile.TemporaryDirectory() as d:
            tsv = os.path.join(d, "s.tsv")
            with open(tsv, "w") as f:
                f.write("species\tgene\tsymbol\nhsa\tENSG00000234456\tMAGI2-AS3\nmmu\tENSMUSG1\tMalat1\n")
            db = os.path.join(d, "ref.db")
            conn = sqlite3.connect(db)
            conn.execute("CREATE TABLE lncrna_reference (species TEXT, transcript_id TEXT, gene TEXT)")
            conn.execute("INSERT INTO lncrna_reference VALUES ('hsa_hg38', 'ENST00000616796', 'ENSG00000234456')")
            conn.commit()
            conn.close()
            s = LncrnaSymbols("hsa_hg38", db_path=db, tsv_path=tsv)
            self.assertEqual(s.symbol_for("ENST00000616796.4"), "MAGI2-AS3")
            self.assertEqual(s.symbol_for("ENSG00000234456"), "MAGI2-AS3")
            self.assertIsNone(s.symbol_for("ENST00000000000"))
            self.assertIsNone(LncrnaSymbols("mmu", db_path=db, tsv_path=tsv).symbol_for("ENSG00000234456"))
            self.assertEqual(LncrnaSymbols("hsa_hg19", db_path=db, tsv_path=tsv).symbol_for("ENSG00000234456"), "MAGI2-AS3")


if __name__ == "__main__":
    unittest.main()
