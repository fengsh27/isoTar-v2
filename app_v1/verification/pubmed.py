"""PubMed E-utilities: query terms, a rate-limited client and a shared SQLite cache.

The query design was validated on 221,035 LUAD/LUSC triplets (isotar-justification
repository). The cache key is the public query string, so any change to a term
below invalidates the cache for that arm and nothing else.

Per arm (mir_gene, mir_lnc, gene_lnc) three nested queries on title/abstract:
    tier 1: A AND B AND <tissue cancer terms>      (not run when tissue is "none")
    tier 2: A AND B AND <any cancer terms>
    tier 3: A AND B
Per triplet with a symbol-bearing lncRNA, a same-paper query: miRNA AND gene AND lncRNA.
Editorial notices (expression of concern, retraction notice, retracted publication)
are excluded by publication type.
"""
import http.client
import json
import os
import re
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL_NAME = "isotar-v2-verification"
DEFAULT_RETMAX = 20

TISSUE = {
    "lung": ('("lung cancer"[tiab] OR "lung neoplasms"[mh] OR "non-small cell lung"[tiab] OR "non-small-cell lung"[tiab] '
             'OR NSCLC[tiab] OR "lung adenocarcinoma"[tiab] OR "lung squamous"[tiab] OR LUAD[tiab] OR LUSC[tiab] '
             'OR "lung carcinoma"[tiab] OR "pulmonary adenocarcinoma"[tiab] OR "lung tumor"[tiab] OR "lung tumour"[tiab])'),
    "breast": '("breast cancer"[tiab] OR "breast neoplasms"[mh] OR "breast carcinoma"[tiab] OR "breast tumor"[tiab] OR "breast tumour"[tiab])',
    "colorectal": '("colorectal cancer"[tiab] OR "colorectal neoplasms"[mh] OR "colon cancer"[tiab] OR "rectal cancer"[tiab] OR CRC[tiab])',
    "gastric": '("gastric cancer"[tiab] OR "stomach neoplasms"[mh] OR "stomach cancer"[tiab] OR "gastric carcinoma"[tiab])',
    "liver": '("hepatocellular carcinoma"[tiab] OR "liver neoplasms"[mh] OR "liver cancer"[tiab] OR HCC[tiab])',
    "prostate": '("prostate cancer"[tiab] OR "prostatic neoplasms"[mh] OR "prostate carcinoma"[tiab])',
    "pancreatic": '("pancreatic cancer"[tiab] OR "pancreatic neoplasms"[mh] OR "pancreatic ductal adenocarcinoma"[tiab] OR PDAC[tiab])',
    "glioma": '(glioma*[tiab] OR glioblastoma*[tiab] OR "brain neoplasms"[mh] OR GBM[tiab])',
}
CANCER = ('(cancer*[tiab] OR tumor*[tiab] OR tumour*[tiab] OR carcinoma*[tiab] OR neoplasm*[tiab] OR oncogen*[tiab] '
          'OR malignan*[tiab] OR glioma*[tiab] OR leukemia*[tiab] OR lymphoma*[tiab] OR melanoma*[tiab] '
          'OR sarcoma*[tiab] OR neoplasms[mh])')
NOTICES = '("expression of concern"[pt] OR "retraction notice"[pt] OR "retracted publication"[pt])'

TIERS = ("1", "2", "3")
ARM_KINDS = ("mir_gene", "mir_lnc", "gene_lnc")


class PubMedError(Exception):
    """E-utilities could not be reached after retries. The message never
    contains the API key -- only the endpoint and the public query."""


def mirna_family(name):
    """hsa-miR-141-5p -> miR-141 ; hsa-let-7c-5p -> let-7c ; mmu-miR-17 -> miR-17.

    Variant headers ("hsa-miR-21-5p,8:A>U,modified") are reduced to the base id
    first: the literature cannot distinguish an in-silico isomiR from its WT."""
    n = name.strip().split(",")[0]
    n = re.sub(r"^[a-z]{3}-", "", n)
    n = re.sub(r"-(3p|5p)$", "", n)
    return n


def mirna_term(name):
    """Family-level disjunction of the common notations for a miRNA."""
    fam = mirna_family(name)
    species = re.match(r"^([a-z]{3})-", name.strip())
    sp = species.group(1) if species else "hsa"
    if fam.lower().startswith("mir-"):
        num = fam[4:]
        forms = ["miR-" + num, "microRNA-" + num, "miRNA-" + num, sp + "-miR-" + num,
                 "miR " + num, "microRNA " + num]
    else:  # let-7 family and anything else
        forms = [fam, sp + "-" + fam]
    return "(" + " OR ".join('"%s"[tiab]' % f for f in forms) + ")"


def sym_term(sym):
    return '"%s"[tiab]' % sym


def is_unsearchable(symbol):
    """An lncRNA with no HGNC symbol (bare Ensembl id, or nothing) cannot be
    searched by name: Ensembl ids never appear in abstracts."""
    s = (symbol or "").strip()
    return (not s) or s.startswith("ENSG") or s.startswith("ENST") or s.startswith("ENSMUS")


def resolve_tissue(value):
    """Turn the request's tissue parameter into (terms, label).

    Accepts a preset name (see TISSUE), "none" (tier 1 not run), or a raw
    PubMed boolean containing [tiab]/[mh] fields. Raises ValueError otherwise."""
    v = (value or "lung").strip()
    low = v.lower()
    if low in TISSUE:
        return TISSUE[low], low
    if low == "none":
        return None, "none"
    if "[" in v and "]" in v:
        return v, "custom"
    raise ValueError("tissue must be one of %s, 'none', or a raw PubMed boolean with [tiab]/[mh] fields"
                     % ", ".join(sorted(TISSUE)))


def arm_entities(kind, mirna, gene, lncrna):
    """The two search terms of an arm."""
    if kind == "mir_gene":
        return mirna_term(mirna), sym_term(gene)
    if kind == "mir_lnc":
        return mirna_term(mirna), sym_term(lncrna)
    return sym_term(gene), sym_term(lncrna)


def tier_queries(x, y, tissue_terms, exclude_notices=True):
    """[(tier, query or None)] -- None when the tier is not run."""
    base = "%s AND %s" % (x, y)
    suffix = " NOT %s" % NOTICES if exclude_notices else ""
    out = []
    for tier, extra in (("1", tissue_terms), ("2", CANCER), ("3", None)):
        if tier == "1" and extra is None:
            out.append((tier, None))
            continue
        q = ("%s AND %s" % (base, extra)) if extra else base
        out.append((tier, q + suffix))
    return out


def same_paper_query(mirna, gene, lncrna, exclude_notices=True):
    q = "%s AND %s AND %s" % (mirna_term(mirna), sym_term(gene), sym_term(lncrna))
    return q + (" NOT %s" % NOTICES if exclude_notices else "")


def pubmed_url(pmid):
    return "https://pubmed.ncbi.nlm.nih.gov/%s/" % pmid


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

class PubMedCache(object):
    """Shared SQLite cache of E-utilities responses, keyed by the public query.

    One file serves every job, so an arm shared by many jobs is asked once.
    max_age_days=0 keeps entries forever; otherwise a stale entry is refetched.
    A None path gives an in-memory cache (tests)."""

    def __init__(self, path=None, max_age_days=0):
        self.path = path
        self.max_age = int(max_age_days or 0) * 86400
        if path:
            os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
            self.conn = sqlite3.connect(path, timeout=60)
        else:
            self.conn = sqlite3.connect(":memory:")
        self.conn.execute("CREATE TABLE IF NOT EXISTS cache "
                          "(key TEXT PRIMARY KEY, value TEXT NOT NULL, fetched_at INTEGER NOT NULL)")
        self.conn.commit()

    def get(self, key):
        row = self.conn.execute("SELECT value, fetched_at FROM cache WHERE key = ?", (key,)).fetchone()
        if row is None:
            return None
        if self.max_age and time.time() - row[1] > self.max_age:
            return None
        return json.loads(row[0])

    def put(self, key, value):
        self.conn.execute("INSERT OR REPLACE INTO cache (key, value, fetched_at) VALUES (?, ?, ?)",
                          (key, json.dumps(value), int(time.time())))
        self.conn.commit()

    def close(self):
        self.conn.close()


# ---------------------------------------------------------------------------
# Client
# ---------------------------------------------------------------------------

def _default_fetch(url, timeout=60):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read().decode("utf-8")


class PubMedClient(object):
    """Rate-limited E-utilities client.

    Keeps at least `min_interval` seconds between live calls (NCBI allows ten
    per second with an API key). Retries transient failures with a capped
    backoff; a web job should fail fast rather than wait out an outage, so the
    default budget is about a minute. `fetch(url) -> str` is injectable for tests."""

    def __init__(self, api_key, email="", cache=None, fetch=None,
                 min_interval=0.11, retries=6, max_sleep=30):
        self.api_key = api_key
        self.email = email or ""
        self.cache = cache if cache is not None else PubMedCache()
        self.fetch = fetch or _default_fetch
        self.min_interval = min_interval
        self.retries = retries
        self.max_sleep = max_sleep
        self.last = 0.0
        self.live_calls = 0

    def _get(self, endpoint, params):
        public = dict(params)
        key = endpoint + "?" + urllib.parse.urlencode(sorted(public.items()))
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        full = dict(params, api_key=self.api_key, email=self.email, tool=TOOL_NAME)
        url = EUTILS + endpoint + "?" + urllib.parse.urlencode(full)
        err = None
        for attempt in range(self.retries):
            wait = self.min_interval - (time.time() - self.last)
            if wait > 0:
                time.sleep(wait)
            self.last = time.time()
            try:
                data = json.loads(self.fetch(url))
                if "error" in data and "esearchresult" not in data and "result" not in data:
                    raise urllib.error.URLError(str(data["error"]))
                self.cache.put(key, data)
                self.live_calls += 1
                return data
            except (urllib.error.URLError, ValueError, OSError, http.client.HTTPException) as e:
                err = e
                time.sleep(min(2 ** attempt, self.max_sleep))
        # Never include `url` here: it carries the API key.
        raise PubMedError("E-utilities failed after %d attempts: %s [%s %s]"
                          % (self.retries, err, endpoint, str(public.get("term", ""))[:200]))

    def esearch(self, term, retmax=DEFAULT_RETMAX):
        d = self._get("esearch.fcgi", {"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json"})
        r = d.get("esearchresult", {})
        return int(r.get("count", 0) or 0), list(r.get("idlist", []))

    def esummary(self, pmids):
        """{pmid: (title, year, journal)} for the given PMIDs, in chunks of 200."""
        out = {}
        pmids = sorted(set(p for p in pmids if p))
        for i in range(0, len(pmids), 200):
            chunk = pmids[i:i + 200]
            d = self._get("esummary.fcgi", {"db": "pubmed", "id": ",".join(chunk), "retmode": "json"})
            result = d.get("result", {})
            for pid in chunk:
                rec = result.get(pid, {}) or {}
                out[pid] = (rec.get("title", ""), (rec.get("pubdate", "") or "")[:4], rec.get("source", ""))
        return out


# ---------------------------------------------------------------------------
# Search driver (pure apart from the client)
# ---------------------------------------------------------------------------

def distinct_arms(triplets):
    """Ordered set of (kind, miRNA, gene, lncRNA_symbol) arms to search.

    The two lncRNA arms are skipped when the lncRNA has no symbol. Triplets are
    dicts with keys miRNA, gene, lncRNA (the search symbol)."""
    seen = {}
    for t in triplets:
        m, g, l = t["miRNA"], t["gene"], t["lncRNA"]
        seen.setdefault(("mir_gene", m, g, ""), None)
        if not is_unsearchable(l):
            seen.setdefault(("mir_lnc", m, "", l), None)
            seen.setdefault(("gene_lnc", "", g, l), None)
    return list(seen.keys())


def search_arms(client, triplets, tissue_terms, retmax=DEFAULT_RETMAX, exclude_notices=True, progress=None):
    """Run the tiered search for every distinct arm.

    Returns {(kind, m, g, l): {"1": (count or "not_run", [pmids]), "2": ..., "3": ...}}.
    `progress(done, total)` is called every 25 arms and at the end."""
    arms = distinct_arms(triplets)
    results = {}
    for n, arm in enumerate(arms, 1):
        kind, m, g, l = arm
        x, y = arm_entities(kind, m, g, l)
        res = {}
        for tier, query in tier_queries(x, y, tissue_terms, exclude_notices):
            if query is None:
                res[tier] = ("not_run", [])
                continue
            res[tier] = client.esearch(query, retmax)
        results[arm] = res
        if progress and (n % 25 == 0 or n == len(arms)):
            progress(n, len(arms))
    if progress and not arms:
        progress(0, 0)
    return results


def search_same_paper(client, triplets, retmax=DEFAULT_RETMAX, exclude_notices=True):
    """{triplet id: (hit 0/1 or None when unsearchable, count, [pmids])}."""
    out = {}
    for t in triplets:
        m, g, l = t["miRNA"], t["gene"], t["lncRNA"]
        if is_unsearchable(l):
            out[t["id"]] = (None, 0, [])
            continue
        count, ids = client.esearch(same_paper_query(m, g, l, exclude_notices), retmax)
        out[t["id"]] = (1 if count else 0, count, ids)
    return out


def representative_note(count, pmid, titles):
    """'n=<count>; <year> <journal>: <title>' for the first returned record."""
    title, year, journal = titles.get(pmid, ("", "", ""))
    return "n=%s; %s %s: %s" % (count, year, journal, title)
