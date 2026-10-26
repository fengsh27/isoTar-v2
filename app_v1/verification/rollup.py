"""Roll per-arm PubMed tiers up to per-triplet rows. Pure functions.

Per arm the tier is the lowest-numbered tier with a hit (1 tissue cancer, 2 any
cancer, 3 symbols only), "none" when no tier hit, "not_run" when tier 1 was
skipped and nothing else hit, and "unsearchable" for the two lncRNA arms of an
lncRNA that has no symbol. The triplet tier is tier{worst arm} when all three
arms hit, "partial" when one or two hit, "none" otherwise.

Column order matches the published verification tables, then the site-specific
extras (ids and scores from the network job).
"""
from app_v1.verification.pubmed import is_unsearchable, pubmed_url, representative_note

PUBLISHED_COLUMNS = [
    "id", "miRNA", "gene", "lncRNA", "triplet_tier", "arms_supported", "same_paper_hit", "same_paper_url",
    "mir_gene_pair", "mir_gene_tier", "mir_gene_url", "mir_gene_note",
    "mir_lnc_pair", "mir_lnc_tier", "mir_lnc_url", "mir_lnc_note",
    "gene_lnc_pair", "gene_lnc_tier", "gene_lnc_url", "gene_lnc_note",
]
EXTRA_COLUMNS = ["gene_id", "lncrna_id", "mirna_node", "pair_score", "bridge_score", "bridge_is_variant"]
COLUMNS = PUBLISHED_COLUMNS + EXTRA_COLUMNS

SUPPORTED = ("1", "2", "3")


def best_tier(res):
    """(tier, [pmids], count) from {"1": (count|"not_run", ids), ...}."""
    for tier in SUPPORTED:
        count, ids = res.get(tier, (0, []))
        if count != "not_run" and count:
            return tier, ids, count
    if any(res.get(t, (0, []))[0] == "not_run" for t in SUPPORTED):
        return "not_run", [], 0
    return "none", [], 0


def short_mirna(m):
    return m[4:] if m.startswith("hsa-") else m


def rollup_rows(triplets, arm_results, same_results, titles):
    """Build one row per triplet.

    triplets:    dicts with id, miRNA, gene, lncRNA (+ the EXTRA_COLUMNS keys)
    arm_results: output of pubmed.search_arms
    same_results: output of pubmed.search_same_paper
    titles:      {pmid: (title, year, journal)}"""
    rows = []
    for t in triplets:
        m, g, l = t["miRNA"], t["gene"], t["lncRNA"]
        row = {"id": t["id"], "miRNA": m, "gene": g, "lncRNA": l}
        for k in EXTRA_COLUMNS:
            row[k] = t.get(k, "")
        tiers = []
        for kind, key in (("mir_gene", ("mir_gene", m, g, "")),
                          ("mir_lnc", ("mir_lnc", m, "", l)),
                          ("gene_lnc", ("gene_lnc", "", g, l))):
            if kind != "mir_gene" and is_unsearchable(l):
                tier, url, note = "unsearchable", "", ""
            else:
                res = arm_results.get(key)
                if res is None:
                    tier, url, note = "", "", ""
                else:
                    tier, ids, count = best_tier(res)
                    if tier in SUPPORTED and ids:
                        url, note = pubmed_url(ids[0]), representative_note(count, ids[0], titles)
                    else:
                        url, note = "", ""
            row[kind + "_tier"], row[kind + "_url"], row[kind + "_note"] = tier, url, note
            tiers.append(tier)
        s = short_mirna(m)
        row["mir_gene_pair"] = "%s to %s" % (s, g)
        row["mir_lnc_pair"] = "%s to %s" % (s, l)
        row["gene_lnc_pair"] = "%s to %s" % (g, l)
        supported = [x for x in tiers if x in SUPPORTED]
        n = len(supported)
        row["triplet_tier"] = ("tier%s" % max(supported)) if n == 3 else ("none" if n == 0 else "partial")
        row["arms_supported"] = n
        hit, count, ids = same_results.get(t["id"], (None, 0, []))
        if hit is None:
            row["same_paper_hit"], row["same_paper_url"] = "", ""
        else:
            row["same_paper_hit"] = hit
            row["same_paper_url"] = pubmed_url(ids[0]) if hit and ids else ""
        rows.append(row)
    rows.sort(key=lambda r: int(r["id"]))
    return rows


def summarize(rows):
    """Counts the UI shows: triplets, per arms_supported, per tier, same-paper hits."""
    out = {
        "triplets": len(rows),
        "with_any_arm": sum(1 for r in rows if r["arms_supported"] >= 1),
        "with_two_arms": sum(1 for r in rows if r["arms_supported"] == 2),
        "with_three_arms": sum(1 for r in rows if r["arms_supported"] == 3),
        "same_paper_hits": sum(1 for r in rows if r["same_paper_hit"] == 1),
        "unsearchable_lncrna": sum(1 for r in rows if r["mir_lnc_tier"] == "unsearchable"),
        "triplet_tier": {},
    }
    for r in rows:
        out["triplet_tier"][r["triplet_tier"]] = out["triplet_tier"].get(r["triplet_tier"], 0) + 1
    return out
