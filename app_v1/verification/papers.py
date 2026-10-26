"""Fetch a paper's abstract and, when open access, its full text; pick the
sentences that mention given terms. Standard library; runs under 3.6 and 3.11.

Full text comes from PMC via E-utilities efetch; if the publisher blocks XML
download there, Europe PMC's fullTextXML endpoint is tried.
"""
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
EUROPE_PMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/"


def _get(url, tries=5, timeout=60):
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read().decode("utf-8")
        except (urllib.error.HTTPError, urllib.error.URLError, OSError) as e:
            if attempt == tries - 1:
                raise
            time.sleep(1 + 2 * attempt)


def _text(el):
    return " ".join(el.itertext()) if el is not None else ""


def fetch_paper(pmid, api_key="", email="", tool="isotar-v2-verification", fetch=None):
    """{pmid, title, journal, year, pmcid, abstract, body, source, found}.

    `fetch(url) -> str` is injectable for tests. The API key goes into the URL
    only; callers must not log the URL."""
    get = fetch or _get
    auth = "&api_key=%s&email=%s&tool=%s" % (urllib.parse.quote(api_key or ""), urllib.parse.quote(email or ""), tool)
    xml = get("%sefetch.fcgi?db=pubmed&id=%s&retmode=xml%s" % (EUTILS, pmid, auth))
    art = ET.fromstring(xml)
    a = art.find(".//PubmedArticle")
    out = {"pmid": str(pmid), "title": "", "journal": "", "year": "", "pmcid": None,
           "abstract": "", "body": "", "source": "", "found": a is not None}
    if a is None:
        return out
    out["title"] = _text(a.find(".//ArticleTitle")).strip()
    out["journal"] = _text(a.find(".//Journal/ISOAbbreviation")).strip()
    out["year"] = (_text(a.find(".//JournalIssue/PubDate/Year")) or _text(a.find(".//JournalIssue/PubDate/MedlineDate"))).strip()
    out["abstract"] = " ".join(_text(x) for x in a.findall(".//Abstract/AbstractText")).strip()
    out["pmcid"] = next((i.text for i in a.findall(".//ArticleId") if i.get("IdType") == "pmc"), None)
    pmc = out["pmcid"]
    if pmc:
        try:
            xml = get("%sefetch.fcgi?db=pmc&id=%s&retmode=xml%s" % (EUTILS, pmc, auth))
            if "does not allow downloading" not in xml:
                out["body"] = _text(ET.fromstring(xml).find(".//body"))
                out["source"] = "PMC efetch"
        except Exception:  # noqa: BLE001 -- full text is optional
            pass
        if not out["body"]:
            try:
                xml = get("%s%s/fullTextXML" % (EUROPE_PMC, pmc))
                out["body"] = _text(ET.fromstring(xml).find(".//body"))
                out["source"] = "Europe PMC"
            except Exception:  # noqa: BLE001
                pass
    return out


def split_sentences(text):
    return [re.sub(r"\s+", " ", s).strip() for s in re.split(r"(?<=[.!?])\s+", text or "") if s.strip()]


def sentences_with(paper, terms, per_term=15, max_len=500):
    """{term: [sentences]} from abstract + body, case-insensitive substring match."""
    sents = split_sentences(paper.get("abstract", "") + " " + paper.get("body", ""))
    out = {}
    for t in terms:
        low = t.lower()
        hits = [s[:max_len] for s in sents if low in s.lower()]
        out[t] = hits[:per_term]
    return out
