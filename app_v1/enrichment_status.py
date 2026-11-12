"""State of a job's enrichment run, kept beside its outputs in the job dir.

Enrichment runs as a Celery task (queue ISOTAR_ENRICHMENT_QUEUE), so the API
answers POST at once and the client polls GET. The run's state lives in
<job>/enrichment_status.json:

    {"status": "queued" | "running" | "succeeded" | "failed",
     "params": {"organism", "cutoff", "gene_count"},
     "queued_at", "started_at", "finished_at", "updated_at",
     "error", "stderr"}

Jobs enriched before this file existed have gene_list.csv but no status; they
read as succeeded so their results keep showing. A queued/running state that
has not been updated for STALE_AFTER seconds reads as failed: the worker that
owned it is gone (restart, crash), and leaving it "running" would block a
re-run forever.
"""
import glob
import json
import os
import time

STATUS_FILENAME = "enrichment_status.json"
ACTIVE = ("queued", "running")
STALE_AFTER = int(os.environ.get("ISOTAR_ENRICHMENT_STALE_SECONDS", "3600"))

# Outputs of v2/enrichment_analysis.py in the job dir, plus the cached result
# zip that bundles them.
_ARTIFACT_PATTERNS = ("*.enrichr.reports.txt", "enrichment_dotplot.png", "result.zip")


def status_path(job_dir):
    return os.path.join(job_dir, STATUS_FILENAME)


def write_status(job_dir, data):
    """Atomically replace the status file, stamping updated_at."""
    data["updated_at"] = int(time.time())
    tmp = status_path(job_dir) + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
    os.replace(tmp, status_path(job_dir))
    return data


def update_status(job_dir, **changes):
    """Merge changes into the current status and write it back."""
    data = _read_file(job_dir) or {}
    data.update(changes)
    return write_status(job_dir, data)


def _read_file(job_dir):
    path = status_path(job_dir)
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r") as f:
            return json.load(f)
    except (ValueError, IOError):
        return None


def read_status(job_dir, now=None):
    """The run's state as the API reports it, or None when enrichment has
    never been run for this job."""
    data = _read_file(job_dir)
    if data is None:
        if os.path.exists(os.path.join(job_dir, "gene_list.csv")):
            return {"status": "succeeded", "legacy": True}
        return None
    if is_stale(data, now):
        data = dict(data)
        data["status"] = "failed"
        data["error"] = (
            "enrichment did not finish within {} minutes; the worker running it "
            "may have restarted. Run it again.".format(STALE_AFTER // 60))
    return data


def is_active(data):
    return bool(data) and data.get("status") in ACTIVE


def is_stale(data, now=None):
    if not is_active(data):
        return False
    now = time.time() if now is None else now
    last = data.get("updated_at") or data.get("queued_at") or 0
    return now - last > STALE_AFTER


def clear_artifacts(job_dir):
    """Remove the previous run's outputs so a new run never shows, or zips,
    results it did not produce. Returns the removed paths."""
    removed = []
    for pattern in _ARTIFACT_PATTERNS:
        for path in glob.glob(os.path.join(job_dir, pattern)):
            os.remove(path)
            removed.append(path)
    return removed
