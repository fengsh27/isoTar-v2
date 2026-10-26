"""Driver for one job's literature verification. Runs inside the Celery worker.

Layout inside the job directory:

    verification/
        progress.json        status, step and counters, polled by the API
        triplets.tsv         the triplets searched (one per pair x bridge miRNA)
        verification.tsv     one row per triplet, published columns first
        verification.json    everything the API serves: params, summary, rows, review
        review_queue.json    hits handed to the model
        verdicts.tsv         the model's verdicts, one per (triplet, scope, PMID)
        llm_review.log       stdout/stderr of the python3.11 review script

Steps 1-5 (triplets, arm search, same-paper, roll-up, queue) are standard
library; step 6 shells out to llm_review.py under python3.11 and is skipped when
the queue is empty or the request left llm_review off.
"""
import csv
import json
import os
import subprocess
import time

from app_v1.verification import pubmed
from app_v1.verification.review import VERDICT_COLUMNS, build_review_queue
from app_v1.verification.rollup import COLUMNS, rollup_rows, summarize
from app_v1.verification.symbols import LncrnaSymbols, species_for_genome

VERIFICATION_DIRNAME = "verification"
QUERY_VERSION = 1  # bump when the PubMed terms or tiers change meaning

_HERE = os.path.dirname(os.path.abspath(__file__))
LLM_SCRIPT = os.path.join(_HERE, "llm_review.py")


def verification_dir(job_dir):
    return os.path.join(job_dir, VERIFICATION_DIRNAME)


def _path(job_dir, name):
    return os.path.join(verification_dir(job_dir), name)


def read_progress(job_dir):
    path = _path(job_dir, "progress.json")
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r") as f:
            return json.load(f)
    except (ValueError, IOError):
        return None


def write_progress(job_dir, data):
    os.makedirs(verification_dir(job_dir), exist_ok=True)
    data["updated_at"] = int(time.time())
    tmp = _path(job_dir, "progress.json.tmp")
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
    os.replace(tmp, _path(job_dir, "progress.json"))


def read_result(job_dir):
    path = _path(job_dir, "verification.json")
    if not os.path.exists(path):
        return None
    with open(path, "r") as f:
        return json.load(f)


def result_tsv_path(job_dir):
    return _path(job_dir, "verification.tsv")


# ---------------------------------------------------------------------------
# Triplets from the network
# ---------------------------------------------------------------------------

def _base_mirna(node_id):
    return node_id.split(",")[0]


def build_triplets(network, symbols):
    """One triplet per (pair, bridge miRNA family) from a build_network() result.

    Pairs come highest score first; within a pair bridges are taken by score.
    Variants of the same miRNA collapse onto the base id (the literature cannot
    tell an in-silico isomiR from its WT), keeping the best-scoring bridge. The
    id is the rank. `symbols.symbol_for(lncrna_id)` supplies the search symbol;
    an lncRNA without one keeps its Ensembl id and is unsearchable on two arms."""
    out = []
    seen = set()
    for pair in network.get("pairs", []):
        gene_id, lncrna_id = pair["gene"], pair["lncrna"]
        gene = pair.get("gene_label") or gene_id
        lnc = symbols.symbol_for(lncrna_id) if symbols is not None else None
        bridges = sorted(pair.get("bridges", []), key=lambda b: -(b.get("score") or 0.0))
        for b in bridges:
            base = _base_mirna(b["mirna"])
            key = (base, gene_id, lncrna_id)
            if key in seen:
                continue
            seen.add(key)
            out.append({
                "id": len(out) + 1,
                "miRNA": base,
                "gene": gene,
                "lncRNA": lnc or lncrna_id,
                "gene_id": gene_id,
                "lncrna_id": lncrna_id,
                "mirna_node": b["mirna"],
                "pair_score": "" if pair.get("score") is None else pair["score"],
                "bridge_score": b.get("score", ""),
                "bridge_is_variant": 1 if b.get("is_variant") else 0,
            })
    return out


# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------

def _write_tsv(path, columns, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _read_verdicts(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return [dict(r) for r in csv.DictReader(f, delimiter="\t")]


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def run_verification(job_dir, network, genome, params, config, task_id=None, log=None):
    """Run the whole verification for one job and write its files.

    params: {"tissue": str, "llm_review": bool}
    config: {"api_key", "email", "cache_path", "cache_max_age_days", "retmax",
             "python", "model", "max_papers", "verdict_cache", "llm_timeout",
             "symbols_tsv", "reference_db"}
    Raises on a failure of the search stage after recording it in progress.json;
    a failure of the model stage is recorded in the result and does not raise."""
    log = log or (lambda *a: None)
    started = int(time.time())
    tissue_terms, tissue_label = pubmed.resolve_tissue(params.get("tissue"))
    llm_wanted = bool(params.get("llm_review"))
    progress = {"status": "running", "step": "triplets", "task_id": task_id, "started_at": started,
                "params": {"tissue": tissue_label, "llm_review": llm_wanted},
                "arms_done": 0, "arms_total": 0, "triplets": 0, "queue_total": 0, "reviewed": 0,
                "error": None}
    write_progress(job_dir, progress)

    try:
        species = species_for_genome(genome)
        symbols = LncrnaSymbols(species, db_path=config.get("reference_db"), tsv_path=config.get("symbols_tsv"))
        triplets = build_triplets(network, symbols)
        # A discovery-mode network can hold hundreds of pairs; cap the search
        # at the top-ranked triplets so one job cannot monopolise the worker.
        max_triplets = int(config.get("max_triplets") or 0)
        truncated = bool(max_triplets and len(triplets) > max_triplets)
        if truncated:
            triplets = triplets[:max_triplets]
        progress["triplets"] = len(triplets)
        _write_tsv(_path(job_dir, "triplets.tsv"),
                   ["id", "miRNA", "gene", "lncRNA", "gene_id", "lncrna_id", "mirna_node", "pair_score",
                    "bridge_score", "bridge_is_variant"], triplets)

        cache = pubmed.PubMedCache(config.get("cache_path"), config.get("cache_max_age_days", 0))
        client = pubmed.PubMedClient(config["api_key"], config.get("email", ""), cache=cache)
        retmax = int(config.get("retmax") or pubmed.DEFAULT_RETMAX)

        progress["step"] = "search_arms"
        progress["arms_total"] = len(pubmed.distinct_arms(triplets))
        write_progress(job_dir, progress)

        def on_arm(done, total):
            progress["arms_done"], progress["arms_total"] = done, total
            write_progress(job_dir, progress)

        arm_results = pubmed.search_arms(client, triplets, tissue_terms, retmax=retmax, progress=on_arm)

        progress["step"] = "same_paper"
        write_progress(job_dir, progress)
        same_results = pubmed.search_same_paper(client, triplets, retmax=retmax)

        progress["step"] = "summaries"
        write_progress(job_dir, progress)
        pmids = set()
        for res in arm_results.values():
            for tier in pubmed.TIERS:
                count, ids = res.get(tier, (0, []))
                if ids:
                    pmids.add(ids[0])
        for hit, count, ids in same_results.values():
            if hit and ids:
                pmids.add(ids[0])
        titles = client.esummary(sorted(pmids))
        cache.close()
        log("verification search done job_dir=%s triplets=%d arms=%d live_calls=%d",
            job_dir, len(triplets), len(arm_results), client.live_calls)

        progress["step"] = "rollup"
        write_progress(job_dir, progress)
        rows = rollup_rows(triplets, arm_results, same_results, titles)
        _write_tsv(result_tsv_path(job_dir), COLUMNS, rows)
        summary = summarize(rows)
        summary["truncated"] = truncated
        summary["max_triplets"] = max_triplets

        queue = build_review_queue(rows, tissue_tier="1" if tissue_terms else "")
        with open(_path(job_dir, "review_queue.json"), "w") as f:
            json.dump(queue, f, indent=2)
        progress["queue_total"] = len(queue)
    except Exception as e:  # noqa: BLE001 -- recorded, then re-raised for the task
        progress["status"] = "failed"
        progress["step"] = "failed"
        progress["error"] = str(e)
        progress["finished_at"] = int(time.time())
        write_progress(job_dir, progress)
        raise

    review = {"status": "disabled" if not llm_wanted else ("skipped" if not queue else "pending"),
              "model": config.get("model"), "queue_size": len(queue), "reviewed": 0,
              "cached": 0, "errors": 0, "verdicts": [], "error": None}
    if llm_wanted and queue:
        progress["step"] = "review"
        write_progress(job_dir, progress)
        review.update(_run_llm_review(job_dir, config, tissue_label, log))
        progress["reviewed"] = review.get("reviewed", 0)

    result = {
        "query_version": QUERY_VERSION,
        "params": {"tissue": tissue_label, "tissue_terms": tissue_terms, "llm_review": llm_wanted},
        "genome": genome,
        "started_at": started,
        "finished_at": int(time.time()),
        "summary": summary,
        "columns": COLUMNS,
        "rows": rows,
        "queue": queue,
        "review": review,
    }
    tmp = _path(job_dir, "verification.json.tmp")
    with open(tmp, "w") as f:
        json.dump(result, f)
    os.replace(tmp, _path(job_dir, "verification.json"))

    progress["status"] = "succeeded"
    progress["step"] = "done"
    progress["finished_at"] = result["finished_at"]
    write_progress(job_dir, progress)
    return result


def _run_llm_review(job_dir, config, tissue_label, log):
    """Run llm_review.py under python3.11; never raises. Returns the review dict fields."""
    out = {"status": "failed", "verdicts": [], "reviewed": 0, "cached": 0, "errors": 0, "error": None}
    verdicts_path = _path(job_dir, "verdicts.tsv")
    cmd = [config.get("python") or "python3.11", "-u", config.get("llm_script") or LLM_SCRIPT,
           _path(job_dir, "review_queue.json"), verdicts_path,
           "--model", str(config.get("model") or ""),
           "--max-papers", str(config.get("max_papers") or 40),
           "--tissue", tissue_label]
    if config.get("verdict_cache"):
        cmd += ["--cache", config["verdict_cache"]]
    log_path = _path(job_dir, "llm_review.log")
    try:
        with open(log_path, "w") as logf:
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=logf, universal_newlines=True,
                                  timeout=int(config.get("llm_timeout") or 3600))
    except subprocess.TimeoutExpired:
        out["error"] = "model review timed out"
        log("verification review timed out job_dir=%s", job_dir)
        return out
    except OSError as e:
        out["error"] = "could not start model review: %s" % e
        log("verification review failed to start job_dir=%s error=%s", job_dir, e)
        return out

    # The script prints one JSON line last; anything before it is progress.
    stats = {}
    for line in (proc.stdout or "").strip().splitlines()[::-1]:
        line = line.strip()
        if line.startswith("{"):
            try:
                stats = json.loads(line)
            except ValueError:
                stats = {}
            break
    out["verdicts"] = _read_verdicts(verdicts_path)
    for k in ("reviewed", "cached", "errors"):
        out[k] = int(stats.get(k, 0) or 0)
    if proc.returncode != 0:
        tail = ""
        try:
            with open(log_path, "r") as f:
                tail = f.read()[-2000:].strip()
        except IOError:
            pass
        out["error"] = stats.get("error") or ("model review exited with code %d: %s" % (proc.returncode, tail[-400:]))
        log("verification review failed job_dir=%s rc=%s", job_dir, proc.returncode)
        return out
    out["status"] = "succeeded"
    log("verification review done job_dir=%s reviewed=%s cached=%s errors=%s",
        job_dir, out["reviewed"], out["cached"], out["errors"])
    return out
