"""Select the PubMed hits worth reading. Pure functions.

A keyword co-occurrence is not a reported interaction, so only the hits that
could change a conclusion are read by the model: every same-paper hit, every
triplet with three supported arms, and two-arm triplets with a tissue-tier arm
or an experimental-looking representative title.
"""
import re

from app_v1.verification.pubmed import is_unsearchable, mirna_family

EXPERIMENTAL = re.compile(
    r"luciferase|reporter|spong|axis|directly target|target of|pull-?down|immunoprecipitation|\bRIP\b|ceRNA|competing endogenous",
    re.I)
REASON_ORDER = {"same_paper": 0, "three_arms": 1, "two_arms_tissue": 2, "two_arms_experimental": 3}
VERDICT_COLUMNS = ["id", "scope", "pmid", "verdict", "arm_match", "tissue", "evidence_type", "sentences", "note"]
ARMS = ("mir_gene", "mir_lnc", "gene_lnc")


def pmid_of(url):
    m = re.search(r"/(\d+)/?$", url or "")
    return m.group(1) if m else ""


def build_review_queue(rows, tissue_tier="1"):
    """[{id, miRNA, gene, lncRNA, reason, arms_supported, triplet_tier,
         pmids: [[scope, pmid], ...], terms: [..]}] ordered by reason then id."""
    out = []
    for r in rows:
        n = int(r["arms_supported"])
        reasons = []
        if str(r.get("same_paper_hit")) == "1":
            reasons.append("same_paper")
        if n == 3:
            reasons.append("three_arms")
        if n == 2:
            tiers = [r[k + "_tier"] for k in ARMS]
            notes = " ".join(r[k + "_note"] for k in ARMS)
            if tissue_tier in tiers:
                reasons.append("two_arms_tissue")
            elif EXPERIMENTAL.search(notes):
                reasons.append("two_arms_experimental")
        if not reasons:
            continue
        pmids = []
        if str(r.get("same_paper_hit")) == "1":
            p = pmid_of(r.get("same_paper_url"))
            if p:
                pmids.append(["same_paper", p])
        for k in ARMS:
            if r[k + "_tier"] in ("1", "2", "3"):
                p = pmid_of(r[k + "_url"])
                if p:
                    pmids.append([k, p])
        terms = [mirna_family(r["miRNA"]), r["gene"]]
        if not is_unsearchable(r["lncRNA"]):
            terms.append(r["lncRNA"])
        out.append({
            "id": r["id"], "miRNA": r["miRNA"], "gene": r["gene"], "lncRNA": r["lncRNA"],
            "reason": ";".join(reasons), "arms_supported": n, "triplet_tier": r["triplet_tier"],
            "pmids": pmids, "terms": terms,
        })
    out.sort(key=lambda q: (min(REASON_ORDER[x] for x in q["reason"].split(";")), int(q["id"])))
    return out
