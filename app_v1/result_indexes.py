"""Build a finished job's result databases ahead of the first request.

The result route builds `result.db` / `lncrna_result.db` / `network.db` lazily
and caches them in the output directory. Done on the request path, that parse
holds the API's single gunicorn worker for its whole duration and, for a large
job, outlives the frontend proxy's 30 s cap, so the first view of a big job
failed with a gateway page while the build kept running server-side.

`build_result_indexes()` runs the same builders from the Celery worker as the
last step of `run_job`, so every `/result` call is a cached read. The lazy
build in the route is unchanged and remains the fallback for jobs that
finished before this existed and for anything the worker could not build.

The module doubles as the post-deploy warm-up for existing jobs:

    python3.6 -m app_v1.result_indexes --dry-run          # list stale indexes
    python3.6 -m app_v1.result_indexes                    # rebuild them all
    python3.6 -m app_v1.result_indexes --job <id> ...     # only these jobs
    python3.6 -m app_v1.result_indexes --since 30         # finished in 30 days

Run it from the repo root (the directory holding `app_v1/`); inside the
container that is `/`. It is a separate process, so it never holds the API
worker, and each database is written atomically, so a concurrent lazy build
from the route is harmless.
"""
import argparse
import json
import os
import sys
import time

from app_v1.lncrna_results import LNCRNA_DB_FILENAME, ensure_lncrna_db
from app_v1.network_results import NETWORK_DB_FILENAME, ensure_network_db
from app_v1.result_db import DB_FILENAME, SCHEMA_VERSION, ensure_db, schema_version

NETWORK_WORKFLOW = "mir-network"
LNCRNA_WORKFLOW = "mir-lncrna"

STATUS_SUCCEEDED = "succeeded"


def planned_indexes(output_dir, workflow):
    """The indexes a workflow's output needs, in build order.

    Each entry is a dict: `name` (for logs), `directory` (the builder's
    argument), `filename`, `builder` (the ensure_* function, looked up by name
    at call time so tests can stub it) and `versioned` (True when the file
    carries result_db.SCHEMA_VERSION in PRAGMA user_version).

    A network pool that was never run (e.g. a TargetScan-only job has no
    lncRNA pool) is left out rather than reported as a failure.
    """
    if workflow == LNCRNA_WORKFLOW:
        return [_spec("lncrna_result.db", output_dir, LNCRNA_DB_FILENAME,
                      "ensure_lncrna_db", True)]
    if workflow == NETWORK_WORKFLOW:
        plan = []
        gene_dir = os.path.join(output_dir, "gene")
        lncrna_dir = os.path.join(output_dir, "lncrna")
        if os.path.isdir(gene_dir):
            plan.append(_spec("gene/result.db", gene_dir, DB_FILENAME,
                              "ensure_db", True))
        if os.path.isdir(lncrna_dir):
            plan.append(_spec("lncrna/lncrna_result.db", lncrna_dir,
                              LNCRNA_DB_FILENAME, "ensure_lncrna_db", True))
        plan.append(_spec("network.db", output_dir, NETWORK_DB_FILENAME,
                          "ensure_network_db", False))
        return plan
    return [_spec("result.db", output_dir, DB_FILENAME, "ensure_db", True)]


def _spec(name, directory, filename, builder, versioned):
    return {
        "name": name,
        "directory": directory,
        "filename": filename,
        "builder": builder,
        "versioned": versioned,
    }


def index_path(spec):
    return os.path.join(spec["directory"], spec["filename"])


def is_current(spec):
    """True when the index exists and (if versioned) is on the current schema,
    i.e. the ensure_* builder would return without rebuilding it."""
    path = index_path(spec)
    if not os.path.exists(path):
        return False
    if spec["versioned"]:
        return schema_version(path) >= SCHEMA_VERSION
    return True


def build_result_indexes(output_dir, workflow, logger=None, job_id=None):
    """Build every index the workflow needs. Never raises.

    Returns one dict per planned index: `name`, `path`, `seconds`, `rebuilt`
    (False when the file was already current and the builder returned at
    once) and `error` (None, or the exception text). A failing index is
    logged and does not stop the others.
    """
    results = []
    for spec in planned_indexes(output_dir, workflow):
        builder = globals()[spec["builder"]]
        rebuilt = not is_current(spec)
        started = time.time()
        error = None
        try:
            builder(spec["directory"])
        except Exception as e:  # noqa: BLE001 - reported, never fatal
            error = str(e) or e.__class__.__name__
        seconds = time.time() - started
        results.append({
            "name": spec["name"],
            "path": index_path(spec),
            "seconds": seconds,
            "rebuilt": rebuilt and error is None,
            "error": error,
        })
        if logger is not None:
            if error is not None:
                logger.warning(
                    "result index failed job_id=%s index=%s seconds=%.1f error=%s",
                    job_id, spec["name"], seconds, error)
            elif rebuilt:
                logger.info("result index built job_id=%s index=%s seconds=%.1f",
                            job_id, spec["name"], seconds)
            else:
                logger.info("result index current job_id=%s index=%s",
                            job_id, spec["name"])
    return results


# ---------------------------------------------------------------------------
# Warm-up CLI for existing jobs
# ---------------------------------------------------------------------------

def discover_jobs(jobs_dir, only=None, since_days=None):
    """Succeeded jobs under jobs_dir, newest finished first.

    Each entry: `job_id`, `output_dir`, `workflow`, `finished_at`. Jobs with
    no readable job.json, a status other than succeeded, or an output
    directory that no longer exists are skipped. `only` restricts to those
    ids; `since_days` drops jobs that finished earlier than that many days
    ago.
    """
    jobs = []
    wanted = set(only) if only else None
    cutoff = time.time() - since_days * 86400 if since_days is not None else None
    for job_id in sorted(os.listdir(jobs_dir)):
        if wanted is not None and job_id not in wanted:
            continue
        job_dir = os.path.join(jobs_dir, job_id)
        meta_path = os.path.join(job_dir, "job.json")
        if not os.path.isfile(meta_path):
            continue
        try:
            with open(meta_path, "r") as fh:
                meta = json.load(fh)
        except (OSError, ValueError):
            continue
        if meta.get("status") != STATUS_SUCCEEDED:
            continue
        finished_at = meta.get("finished_at") or 0
        if cutoff is not None and finished_at < cutoff:
            continue
        output_dir = meta.get("result_path") or os.path.join(job_dir, "output")
        if not os.path.isdir(output_dir):
            continue
        jobs.append({
            "job_id": job_id,
            "output_dir": output_dir,
            "workflow": meta.get("workflow", "mir-target"),
            "finished_at": finished_at,
        })
    jobs.sort(key=lambda j: j["finished_at"], reverse=True)
    return jobs


def main(argv=None, out=None):
    out = out or sys.stdout
    parser = argparse.ArgumentParser(
        prog="python3.6 -m app_v1.result_indexes",
        description="Build or refresh the cached result databases of finished "
                    "jobs, so the API never has to build one on request.")
    parser.add_argument("--jobs-dir",
                        default=os.environ.get("ISOTAR_JOB_DIR", "/opt/out/jobs"),
                        help="job storage root (default: $ISOTAR_JOB_DIR or /opt/out/jobs)")
    parser.add_argument("--job", action="append", default=None, metavar="ID",
                        help="only this job id (repeatable)")
    parser.add_argument("--since", type=float, default=None, metavar="DAYS",
                        help="only jobs that finished within the last DAYS days")
    parser.add_argument("--dry-run", action="store_true",
                        help="list the indexes that would be rebuilt; build nothing")
    args = parser.parse_args(argv)

    if not os.path.isdir(args.jobs_dir):
        out.write("error: jobs dir not found: {}\n".format(args.jobs_dir))
        return 2

    jobs = discover_jobs(args.jobs_dir, only=args.job, since_days=args.since)
    stale = built = failed = 0
    slowest = (0.0, None)
    for job in jobs:
        for spec in planned_indexes(job["output_dir"], job["workflow"]):
            if is_current(spec):
                continue
            stale += 1
            if args.dry_run:
                out.write("stale    {} {} {}\n".format(
                    job["job_id"], job["workflow"], spec["name"]))
        if args.dry_run:
            continue
        for r in build_result_indexes(job["output_dir"], job["workflow"]):
            if r["error"] is not None:
                failed += 1
                out.write("FAILED   {} {} {} {:.1f}s {}\n".format(
                    job["job_id"], job["workflow"], r["name"], r["seconds"], r["error"]))
            elif r["rebuilt"]:
                built += 1
                out.write("built    {} {} {} {:.1f}s\n".format(
                    job["job_id"], job["workflow"], r["name"], r["seconds"]))
                if r["seconds"] > slowest[0]:
                    slowest = (r["seconds"], "{} {}".format(job["job_id"], r["name"]))
        out.flush()

    if args.dry_run:
        out.write("{} succeeded job(s), {} stale index(es)\n".format(len(jobs), stale))
    else:
        out.write("{} succeeded job(s), {} index(es) built, {} failed\n".format(
            len(jobs), built, failed))
        if slowest[1]:
            out.write("slowest build: {} in {:.1f}s\n".format(slowest[1], slowest[0]))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
