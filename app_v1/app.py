import csv
import glob
import json
import os
import re
import shutil
import subprocess
import time
import uuid
import zipfile

from flask import Flask, jsonify, request, send_file
from flask_cors import CORS

from app_v1.celery_app import celery_app
from app_v1.limits import MAX_CORES_PER_JOB, validate_cores
from app_v1.logger import get_logger
from app_v1.version import get_version
from app_v1.result_db import ensure_db, query_genes
from app_v1.lncrna_results import ensure_lncrna_db
from app_v1.network_results import (
    build_network,
    load_pairs,
    DEFAULT_TOP_GENES,
    DEFAULT_TOP_LNCRNA,
)
from app_v1.target_resolver import (
    genome_to_species as _genome_to_species,
    resolve_targets as _resolve_targets,
    validate_gene_targets as _validate_gene_targets,
)
from app_v1.lncrna_reference import (
    validate_lncrna_targets as _validate_lncrna_targets,
    normalize_lncrna_id as _normalize_lncrna_id,
)

logger = get_logger()

BASE_DIR = os.environ.get("ISOTAR_JOB_DIR", "/opt/out/jobs")
os.makedirs(BASE_DIR, exist_ok=True)

# Frontend "workflow" -> mirna_predicting.py "--target-type" (target sequence pool).
_WORKFLOW_TARGET_TYPE = {
    "mir-target": "gene",
    "mir-lncrna": "lncrna",
}
# Tools that cannot run against a lncRNA target pool: TargetScan ignores the
# target FASTA and reads its own precomputed 3' UTR datasets. (PITA is fine -- it
# scores whatever FASTA is passed to -utr.) Mirrors LNCRNA_INCOMPATIBLE_TOOLS in
# v2/mirna_predicting.py.
LNCRNA_INCOMPATIBLE_TOOLS = {"Targetscan"}

# The Network workflow runs BOTH target pools (gene + lncRNA) over a list of
# miRNAs in a single job, then joins gene-target <-> miRNA <-> lncRNA-target into
# a tripartite graph. It is handled separately from _WORKFLOW_TARGET_TYPE because
# it is not a single-pool run.
NETWORK_WORKFLOW = "mir-network"
# Compute scales as (miRNAs x full reference x 2 pools); cap the list so a single
# job can't fan out into an unbounded run.
MAX_NETWORK_MIRNAS = int(os.environ.get("ISOTAR_MAX_NETWORK_MIRNAS", "20"))

app = Flask(__name__)

_cors_origins = os.environ.get("ISOTAR_CORS_ORIGINS", "*")
if _cors_origins != "*":
    _cors_origins = [o.strip() for o in _cors_origins.split(",") if o.strip()]
CORS(app, resources={r"/api/*": {"origins": _cors_origins}})

_VERSION = get_version()


@app.route("/api/v1/version", methods=["GET"])
def version():
    """Report the running backend version (matches the Docker image / git tag)."""
    return jsonify({"name": "isotar-v2-backend", "version": _VERSION})


def _job_path(job_id):
    return os.path.join(BASE_DIR, job_id)


def _job_meta_path(job_id):
    return os.path.join(_job_path(job_id), "job.json")


def _progress_path(job_id):
    return os.path.join(_job_path(job_id), "output", "progress.json")


# mir-network runs every tool against two target pools (gene + lncRNA), each in
# its own output subdirectory, so each pool writes its OWN progress.json. There
# is no top-level output/progress.json for a network job -- without merging
# these the API would serve no progress and the UI would show every tool stuck
# "pending". These are the pool subdirectories run_job() passes as `-o`.
_NETWORK_POOLS = ("gene", "lncrna")


def _combine_tool_status(entries):
    """Combine one tool's status across the pools it ran in. A tool is only
    'done' when every pool finished it; 'running' if any pool is running it;
    'pending' otherwise. Timing spans the pools (earliest start, latest finish)."""
    statuses = [e.get("status") for e in entries]
    if "running" in statuses:
        status = "running"
    elif statuses and all(s == "done" for s in statuses):
        status = "done"
    else:
        status = "pending"

    starts = [e.get("started_at") for e in entries if e.get("started_at") is not None]
    started_at = min(starts) if starts else None
    if status == "done":
        finishes = [e.get("finished_at") for e in entries if e.get("finished_at") is not None]
        finished_at = max(finishes) if finishes else None
    else:
        finished_at = None
    return {"status": status, "started_at": started_at, "finished_at": finished_at}


def _merge_pool_progress(pool_files):
    """Merge per-pool progress.json files into a single per-tool view keyed by
    bare tool name (what the frontend matches against job['tools']). Returns None
    when no pool file is readable."""
    by_tool = {}
    order = []
    updated_at = 0
    for path in pool_files:
        try:
            with open(path, "r") as f:
                data = json.load(f)
        except (ValueError, IOError):
            continue
        updated_at = max(updated_at, data.get("updated_at") or 0)
        for tool, info in (data.get("tools_status") or {}).items():
            if tool not in by_tool:
                by_tool[tool] = []
                order.append(tool)
            by_tool[tool].append(info)

    if not by_tool:
        return None

    tools_status = {tool: _combine_tool_status(by_tool[tool]) for tool in order}
    completed = sum(1 for t in tools_status if tools_status[t]["status"] == "done")
    current = next((t for t in order if tools_status[t]["status"] == "running"), None)
    return {
        "total_tools": len(tools_status),
        "completed_tools": completed,
        "current_tool": current,
        "tools_status": tools_status,
        "updated_at": updated_at,
    }


def _write_meta(job_id, data):
    with open(_job_meta_path(job_id), "w") as f:
        json.dump(data, f, indent=2)


def _load_meta(job_id):
    with open(_job_meta_path(job_id), "r") as f:
        return json.load(f)


def _load_progress(job_id):
    path = _progress_path(job_id)
    if os.path.exists(path):
        try:
            with open(path, "r") as f:
                return json.load(f)
        except (ValueError, IOError):
            return None

    # No top-level progress.json -- a network job writes one per pool instead.
    output_dir = os.path.join(_job_path(job_id), "output")
    pool_files = [
        os.path.join(output_dir, pool, "progress.json")
        for pool in _NETWORK_POOLS
    ]
    pool_files = [p for p in pool_files if os.path.exists(p)]
    if pool_files:
        return _merge_pool_progress(pool_files)
    return None


# Custom miRNA sequence bounds -- mirror v2/mirna_processing.py
# (MIN_LENGTH / MAX_LENGTH / VALID_NUCLEOTIDES).
_MIRNA_SEQ_MIN_LEN = 17
_MIRNA_SEQ_MAX_LEN = 30
_MIRNA_SEQ_ALPHABET = set("ACGTU")


def _validate_mirna_seq(seq):
    """Validate a user-supplied custom miRNA sequence.

    Returns (normalized_seq, error) where error is (message, status) or None.
    Mirrors mirna_processing.py: 17-30 nt, alphabet A C G T U (case-insensitive).
    """
    if not isinstance(seq, str):
        return None, ("mirna_seq must be a string", 400)
    s = seq.strip().upper()
    if not s:
        return None, ("mirna_seq must not be empty", 400)
    invalid = sorted(set(s) - _MIRNA_SEQ_ALPHABET)
    if invalid:
        return None, ("mirna_seq contains invalid nucleotides {} (allowed: A C G T U)".format(invalid), 400)
    if not (_MIRNA_SEQ_MIN_LEN <= len(s) <= _MIRNA_SEQ_MAX_LEN):
        return None, ("mirna_seq length {} is outside the allowed range ({}-{})".format(
            len(s), _MIRNA_SEQ_MIN_LEN, _MIRNA_SEQ_MAX_LEN), 400)
    return s, None


def _process_single_mirna(job_id, mirna_id, fasta_path, modifications, shift, pre_id,
                          mirna_seq=None):
    """Generate mirna.fa for one miRNA: a user-supplied custom sequence when
    mirna_seq is given, otherwise the miRBase WT (+ optional modifications/shift)."""
    cmd = [
        "python2.7",
        "/opt/v2/mirna_processing.py",
        mirna_id,
        "-o",
        fasta_path,
    ]
    if mirna_seq:
        # Custom miRNA: pass the raw sequence straight through. modifications/
        # shift/pre_id need miRBase precursor context a custom sequence lacks,
        # so they are not applied.
        cmd.extend(["--seq", mirna_seq])
    else:
        for mod in modifications:
            cmd.extend(["-m", mod])
        if shift:
            # Use the --shift=VALUE form (not ["-s", shift]): a negative
            # left-shift like "-7|1" starts with "-" and isn't a valid negative
            # number (the "|1" breaks argparse's negative-number matcher), so
            # space-separated argparse treats it as an unknown option and dies
            # with "expected one argument" (exit 2). The =-joined form binds the
            # value directly and parses cleanly.
            cmd.append("--shift=" + shift)
        if modifications and shift:
            cmd.append("-b")
        if pre_id:
            cmd.extend(["--pre-id", pre_id])
    logger.info("step=processing job_id=%s", job_id)
    subprocess.check_call(cmd)


def _process_mirna_list(job_id, mirna_ids, fasta_path, meta):
    """Process each miRNA (WT only) and concatenate into one mirna.fa.

    Per-miRNA processing failures are recorded in meta['mirna_errors'] and that
    miRNA is skipped; the job only aborts if EVERY miRNA failed. stdin is closed
    so an interactive precursor prompt (multi-precursor miRNA without a pre-id)
    raises EOFError instead of hanging the worker."""
    job_dir = _job_path(job_id)
    parts_dir = os.path.join(job_dir, "mirna_parts")
    os.makedirs(parts_dir, exist_ok=True)
    logger.info("step=processing job_id=%s mirnas=%d", job_id, len(mirna_ids))

    # Per-miRNA precursor selection (from the submit request). A multi-precursor
    # miRNA without an entry here still fails processing and is dropped below.
    pre_ids = meta.get("pre_ids") or {}
    mirna_errors = {}
    part_paths = []
    for idx, mirna_id in enumerate(mirna_ids):
        part = os.path.join(parts_dir, "{}.fa".format(idx))
        cmd = ["python2.7", "/opt/v2/mirna_processing.py", mirna_id, "-o", part]
        pre_id = pre_ids.get(mirna_id)
        if pre_id:
            cmd.extend(["--pre-id", pre_id])
        try:
            with open(os.devnull, "rb") as devnull:
                subprocess.check_call(cmd, stdin=devnull)
            part_paths.append(part)
        except subprocess.CalledProcessError as e:
            mirna_errors[mirna_id] = str(e)
            logger.error("mirna processing failed job_id=%s mirna=%s error=%s",
                         job_id, mirna_id, e)

    if not part_paths:
        raise RuntimeError("all miRNAs failed processing: {}".format(mirna_errors))

    with open(fasta_path, "w") as out:
        for part in part_paths:
            with open(part, "r") as pf:
                data = pf.read()
            out.write(data)
            if data and not data.endswith("\n"):
                out.write("\n")

    if mirna_errors:
        meta["mirna_errors"] = mirna_errors


def _run_prediction_pool(job_id, fasta_path, tools, genome, target_type,
                         pool_output_dir, cores, target_file, tool_errors, label=None):
    """Run the selected tools against one target pool as two subprocess groups:
    the python3.6 tools, then miRmap under python3.11. A failure in one group is
    recorded in tool_errors (keyed '<label>:<tools>' when label is set, or the
    bare tool key otherwise) and never discards the other group's results.
    Returns the number of groups attempted."""
    os.makedirs(pool_output_dir, exist_ok=True)
    groups_run = 0

    other_tools = [t for t in tools if t != "miRmap"]
    if other_tools:
        groups_run += 1
        key = "+".join(other_tools) if label is None else "{}:{}".format(label, "+".join(other_tools))
        logger.info("running python3.6 tools=%s pool=%s job_id=%s", other_tools, label or target_type, job_id)
        cmd = [
            "python3.6",
            "/opt/v2/mirna_predicting.py",
            "-c", str(cores),
            "-i", fasta_path,
            "-t",
        ] + other_tools + ["-g", genome, "-tt", target_type, "-o", pool_output_dir]
        if target_file:
            cmd += ["-tf", target_file]
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError as e:
            tool_errors[key] = str(e)
            logger.error("tool group failed job_id=%s pool=%s tools=%s error=%s",
                         job_id, label or target_type, other_tools, e)

    if "miRmap" in tools:
        groups_run += 1
        key = "miRmap" if label is None else "{}:miRmap".format(label)
        logger.info("running python3.11 tools=['miRmap'] pool=%s job_id=%s", label or target_type, job_id)
        cmd = [
            "python3.11",
            "/opt/v2/mirna_predicting.py",
            "-c", str(cores),
            "-i", fasta_path,
            "-t", "miRmap",
            "-g", genome,
            "-tt", target_type,
            "-o", pool_output_dir,
        ]
        if target_file:
            cmd += ["-tf", target_file]
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError as e:
            tool_errors[key] = str(e)
            logger.error("tool group failed job_id=%s pool=%s tools=['miRmap'] error=%s",
                         job_id, label or target_type, e)

    return groups_run


@celery_app.task(bind=True)
def run_job(self, job_id):
    meta = _load_meta(job_id)
    meta["status"] = "running"
    meta["task_id"] = self.request.id
    meta["started_at"] = int(time.time())
    _write_meta(job_id, meta)

    job_dir = _job_path(job_id)
    output_dir = os.path.join(job_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    try:
        tools = meta["tools"]
        genome = meta.get("genome", "hg38")
        cores = int(meta.get("cores", 1))
        workflow = meta.get("workflow", "mir-target")
        is_network = workflow == NETWORK_WORKFLOW
        fasta_path = os.path.join(job_dir, "mirna.fa")

        # --- processing: generate mirna.fa (one miRNA, or a list for network) ---
        if is_network:
            mirna_ids = meta.get("mirna_ids", [])
            logger.info("network job started job_id=%s mirnas=%d tools=%s genome=%s cores=%s",
                        job_id, len(mirna_ids), tools, genome, cores)
            _process_mirna_list(job_id, mirna_ids, fasta_path, meta)
        else:
            mirna_id = meta["mirna_id"]
            logger.info("job started job_id=%s mirna_id=%s tools=%s genome=%s cores=%s",
                        job_id, mirna_id, tools, genome, cores)
            _process_single_mirna(
                job_id, mirna_id, fasta_path,
                meta.get("modifications", []), meta.get("shift"), meta.get("pre_id"),
                meta.get("mirna_seq"),
            )
        meta["step"] = "processing"
        _write_meta(job_id, meta)

        meta["step"] = "predicting"
        _write_meta(job_id, meta)
        logger.info("step=predicting job_id=%s", job_id)

        # The tools run as two independent subprocess groups per pool (python3.6
        # tools, then miRmap under python3.11). A failure in one group must not
        # discard the results of the others -- each runs in its own try/except
        # (in _run_prediction_pool) and failures are recorded in tool_errors. The
        # job is only marked failed if EVERY group failed (nothing to return);
        # otherwise it succeeds with whatever results were produced and the
        # failed tools noted in meta["tool_errors"].
        tool_errors = {}
        groups_run = 0

        if is_network:
            # Run BOTH target pools over the same miRNA list. TargetScan can only
            # scan the 3' UTR (gene) pool -- drop it on the lncRNA pass. No target
            # filter: the ceRNA pairs filter at result-build time, not here.
            groups_run += _run_prediction_pool(
                job_id, fasta_path, tools, genome, "gene",
                os.path.join(output_dir, "gene"), cores, None, tool_errors, label="gene")
            lncrna_tools = [t for t in tools if t not in LNCRNA_INCOMPATIBLE_TOOLS]
            if lncrna_tools:
                groups_run += _run_prediction_pool(
                    job_id, fasta_path, lncrna_tools, genome, "lncrna",
                    os.path.join(output_dir, "lncrna"), cores, None, tool_errors, label="lncrna")
        else:
            target_type = _WORKFLOW_TARGET_TYPE.get(workflow, "gene")
            groups_run += _run_prediction_pool(
                job_id, fasta_path, tools, genome, target_type,
                output_dir, cores, meta.get("target_file"), tool_errors, label=None)

        if tool_errors and len(tool_errors) >= groups_run:
            # Every tool group failed -- there are no results to return.
            meta["status"] = "failed"
            meta["finished_at"] = int(time.time())
            meta["error"] = "All prediction tools failed: {}".format(tool_errors)
            _write_meta(job_id, meta)
            logger.error("job failed job_id=%s error=%s", job_id, tool_errors)
            return

        meta["status"] = "succeeded"
        meta["finished_at"] = int(time.time())
        meta["result_path"] = output_dir
        if tool_errors:
            meta["tool_errors"] = tool_errors
        _write_meta(job_id, meta)
        logger.info("job succeeded job_id=%s tool_errors=%s", job_id, tool_errors)
    except subprocess.CalledProcessError as e:
        meta["status"] = "failed"
        meta["finished_at"] = int(time.time())
        meta["error"] = "Command failed: {}".format(e)
        _write_meta(job_id, meta)
        logger.error("job failed job_id=%s error=%s", job_id, e)
    except Exception as e:
        meta["status"] = "failed"
        meta["finished_at"] = int(time.time())
        meta["error"] = str(e)
        _write_meta(job_id, meta)
        logger.error("job failed job_id=%s error=%s", job_id, e)


def _submit_network_job(data, tools):
    """Validate, resolve, and enqueue a mir-network job (list of miRNAs + both
    target pools + optional ceRNA (gene, lncRNA) pairs)."""
    mirna_ids = data.get("mirna_ids") or []
    pairs_in = data.get("pairs") or []

    if not tools:
        return jsonify({"error": "tools are required"}), 400
    if (not isinstance(mirna_ids, list) or not mirna_ids
            or not all(isinstance(m, str) and m.strip() for m in mirna_ids)):
        return jsonify({"error": "mirna_ids must be a non-empty list of miRNA ID strings"}), 400

    # Dedupe while preserving the user's order.
    seen = set()
    mirna_ids = [m.strip() for m in mirna_ids
                 if m.strip() not in seen and not seen.add(m.strip())]
    if len(mirna_ids) > MAX_NETWORK_MIRNAS:
        return jsonify({"error": "too many miRNAs for one network job",
                        "max": MAX_NETWORK_MIRNAS, "count": len(mirna_ids)}), 400

    # Optional per-miRNA precursor selection: {"<mirna_id>": "<pre_id>", ...}.
    # Only multi-precursor miRNAs (e.g. hsa-let-7a-5p -> let-7a-1/-2/-3) need an
    # entry; unambiguous ids resolve on their own in mirna_processing.py. Keys
    # must reference a submitted miRNA so a typo surfaces as a 400 instead of
    # being silently ignored (and the miRNA then dropped for missing a pre-id).
    pre_ids_in = data.get("pre_ids") or {}
    if not isinstance(pre_ids_in, dict):
        return jsonify({"error": "pre_ids must be an object mapping miRNA id -> precursor id"}), 400
    pre_ids = {}
    for k, v in pre_ids_in.items():
        if not (isinstance(k, str) and isinstance(v, str) and k.strip() and v.strip()):
            return jsonify({"error": "pre_ids must map non-empty miRNA id strings to non-empty precursor id strings"}), 400
        pre_ids[k.strip()] = v.strip()
    unknown_pre_ids = sorted(k for k in pre_ids if k not in mirna_ids)
    if unknown_pre_ids:
        return jsonify({"error": "pre_ids references miRNAs not in mirna_ids",
                        "unknown": unknown_pre_ids}), 400

    cores, cores_err = validate_cores(data.get("cores"))
    if cores_err:
        msg, status = cores_err
        return jsonify({"error": msg, "max_cores_per_job": MAX_CORES_PER_JOB}), status

    genome = data.get("genome", "hg38")

    if not isinstance(pairs_in, list):
        return jsonify({"error": "pairs must be a list of {gene, lncrna} objects"}), 400

    # Resolve each pair's gene symbol/accession to RefSeq before creating job
    # state, so a 400 leaves no orphan dir. A symbol may expand to several RefSeq.
    resolved_pairs = []
    unresolved_genes = []
    for p in pairs_in:
        if not isinstance(p, dict):
            return jsonify({"error": "each pair must be an object with gene and lncrna"}), 400
        gene = (p.get("gene") or "").strip()
        lncrna = (p.get("lncrna") or "").strip()
        if not gene or not lncrna:
            return jsonify({"error": "each pair requires a non-empty gene and lncrna"}), 400
        refseqs, unresolved = _resolve_targets([gene], genome)
        if unresolved:
            unresolved_genes.extend(unresolved)
            continue
        resolved_pairs.append({
            "gene": gene,
            "gene_refseqs": sorted(refseqs),
            "lncrna": lncrna,
        })
    if unresolved_genes:
        return jsonify({
            "error": "unresolved pair genes",
            "unresolved": unresolved_genes,
            "hint": "gene symbols not found in gene_mapping for species {}".format(
                _genome_to_species(genome)),
        }), 400

    job_id = str(uuid.uuid4())
    job_dir = _job_path(job_id)
    os.makedirs(job_dir, exist_ok=True)

    if resolved_pairs:
        with open(os.path.join(job_dir, "pairs.json"), "w") as f:
            json.dump({"pairs": resolved_pairs}, f, indent=2)

    meta = {
        "job_id": job_id,
        "status": "queued",
        "created_at": int(time.time()),
        "mirna_ids": mirna_ids,
        "tools": tools,
        "workflow": NETWORK_WORKFLOW,
        "genome": genome,
        "cores": cores,
        "pair_count": len(resolved_pairs),
        "pre_ids": pre_ids,
    }
    _write_meta(job_id, meta)
    task = run_job.delay(job_id)
    meta["task_id"] = task.id
    _write_meta(job_id, meta)
    logger.info("network job queued job_id=%s mirnas=%d tools=%s genome=%s pairs=%d",
                job_id, len(mirna_ids), tools, genome, len(resolved_pairs))
    return jsonify({"job_id": job_id, "task_id": task.id}), 202


_VALIDATE_TARGET_TYPES = ("gene", "lncrna")
MAX_VALIDATE_TARGETS = int(os.environ.get("ISOTAR_MAX_VALIDATE_TARGETS", "5000"))


@app.route("/api/v1/targets/validate", methods=["POST"])
def validate_targets():
    """Check whether user-supplied targets exist in our reference data.

    body: {targets:[...], genome:"hg38", target_type:"gene"|"lncrna"}
    ->    {genome, species, target_type, results:[{target,valid,matched_by}],
           valid_count, invalid:[...]}
    'gene' checks HGNC symbols / RefSeq accessions in gene_mapping; 'lncrna'
    checks Ensembl/FlyBase/WormBase lncRNA transcript & gene ids. Read-only:
    creates no job state."""
    data = request.get_json(force=True, silent=True) or {}
    targets = data.get("targets")
    genome = data.get("genome", "hg38")
    target_type = data.get("target_type", "gene")

    if target_type not in _VALIDATE_TARGET_TYPES:
        return jsonify({"error": "target_type must be one of {}".format(
            list(_VALIDATE_TARGET_TYPES))}), 400
    if (not isinstance(targets, list) or not targets
            or not all(isinstance(t, str) for t in targets)):
        return jsonify({"error": "targets must be a non-empty list of strings"}), 400
    if len(targets) > MAX_VALIDATE_TARGETS:
        return jsonify({"error": "too many targets (max {})".format(MAX_VALIDATE_TARGETS),
                        "max_targets": MAX_VALIDATE_TARGETS}), 400

    try:
        if target_type == "lncrna":
            out = _validate_lncrna_targets(targets, genome)
        else:
            out = _validate_gene_targets(targets, genome)
    except Exception:
        logger.exception("target validation failed target_type=%s genome=%s",
                         target_type, genome)
        return jsonify({"error": "internal error during target validation"}), 500

    out["genome"] = genome
    out["target_type"] = target_type
    logger.info("validated %d targets target_type=%s genome=%s valid=%d",
                len(targets), target_type, genome, out["valid_count"])
    return jsonify(out)


@app.route("/api/v1/jobs", methods=["POST"])
def submit_job():
    data = request.get_json(force=True, silent=True) or {}
    tools = data.get("tools", [])
    workflow = data.get("workflow", "mir-target")

    # The Network workflow takes a list of miRNAs (not a single mirna_id) and
    # optional (gene, lncRNA) pairs, so it has its own validation/enqueue path.
    if workflow == NETWORK_WORKFLOW:
        return _submit_network_job(data, tools)

    mirna_id = data.get("mirna_id")
    mirna_seq = data.get("mirna_seq")
    modifications = data.get("modifications", [])
    shift = data.get("shift")
    pre_id = data.get("pre_id")
    targets = data.get("targets")

    if not tools or not mirna_id:
        logger.warning("job rejected: missing tools or mirna_id")
        return jsonify({"error": "tools and mirna_id are required"}), 400

    if workflow not in _WORKFLOW_TARGET_TYPE:
        return jsonify({"error": "workflow must be one of {}".format(
            sorted(_WORKFLOW_TARGET_TYPE))}), 400
    if workflow == "mir-lncrna":
        bad = [t for t in tools if t in LNCRNA_INCOMPATIBLE_TOOLS]
        if bad:
            return jsonify({
                "error": "tools {} are not supported for the miR-LncRNA workflow".format(bad),
                "hint": "TargetScan/PITA only work on 3' UTR (gene) targets",
            }), 400

    if not isinstance(modifications, list):
        logger.warning("job rejected: modifications not a list mirna_id=%s", mirna_id)
        return jsonify({"error": "modifications must be a list of strings"}), 400

    if shift is not None and not isinstance(shift, str):
        logger.warning("job rejected: invalid shift mirna_id=%s shift=%s", mirna_id, shift)
        return jsonify({"error": "shift must be a string in format 'left|right' (e.g. '-4|-6')"}), 400

    # Custom miRNA: a raw sequence supplied alongside the (free-form) mirna_id
    # label. When present it bypasses the miRBase lookup; modifications/shift do
    # not apply to a custom sequence, so they are ignored downstream.
    if mirna_seq is not None:
        # The custom name flows into the FASTA header and downstream output
        # filenames (<header>_<tool>_results.txt), so constrain it to a safe
        # charset -- no spaces or path separators.
        if not isinstance(mirna_id, str) or not re.match(r'^[A-Za-z0-9._-]+$', mirna_id):
            logger.warning("job rejected: invalid custom mirna_id=%r", mirna_id)
            return jsonify({
                "error": "mirna_id must contain only letters, digits, '.', '_' or '-' "
                         "(no spaces or path separators) for a custom miRNA sequence",
            }), 400
        mirna_seq, seq_err = _validate_mirna_seq(mirna_seq)
        if seq_err:
            msg, status = seq_err
            logger.warning("job rejected: %s mirna_id=%s", msg, mirna_id)
            return jsonify({"error": msg}), status

    if targets is not None:
        if not isinstance(targets, list) or not targets or not all(isinstance(g, str) for g in targets):
            return jsonify({"error": "targets must be a non-empty list of target id strings"}), 400

    cores, cores_err = validate_cores(data.get("cores"))
    if cores_err:
        msg, status = cores_err
        logger.warning("job rejected: %s mirna_id=%s cores=%r", msg, mirna_id, data.get("cores"))
        return jsonify({"error": msg, "max_cores_per_job": MAX_CORES_PER_JOB}), status

    genome = data.get("genome", "hg38")
    target_type = _WORKFLOW_TARGET_TYPE[workflow]

    # Resolve/validate targets before creating any job state so a 400 leaves no
    # orphan dir. gene -> RefSeq accessions filtered from the 3' UTR reference;
    # lncrna -> normalized transcript/gene ids filtered from the lncRNA reference.
    # Strict: any unknown target rejects the whole job.
    target_tokens = None
    if targets:
        if target_type == "lncrna":
            vres = _validate_lncrna_targets(targets, genome)
            if vres["invalid"]:
                logger.warning("job rejected: unknown lncRNA targets unknown=%s species=%s",
                               vres["invalid"], vres["species"])
                return jsonify({
                    "error": "unresolved targets",
                    "unresolved": vres["invalid"],
                    "hint": "lncRNA transcript/gene ids not found for species {}".format(
                        vres["species"]),
                }), 400
            target_tokens = sorted({n for n in
                                    (_normalize_lncrna_id(t) for t in targets) if n})
        else:
            accessions, unresolved = _resolve_targets(targets, genome)
            if unresolved:
                logger.warning("job rejected: unresolved targets unresolved=%s species=%s",
                               unresolved, _genome_to_species(genome))
                return jsonify({
                    "error": "unresolved targets",
                    "unresolved": unresolved,
                    "hint": "gene symbols not found in gene_mapping for species {}".format(
                        _genome_to_species(genome)),
                }), 400
            if not accessions:
                return jsonify({"error": "targets resolved to zero RefSeq accessions"}), 400
            target_tokens = sorted(accessions)

    job_id = str(uuid.uuid4())
    job_dir = _job_path(job_id)
    os.makedirs(job_dir, exist_ok=True)

    target_file = None
    if target_tokens:
        target_file = os.path.join(job_dir, "targets.txt")
        with open(target_file, "w") as f:
            for tok in target_tokens:
                f.write(tok + "\n")
        logger.info("resolved %d target tokens to %d %s identifiers job_id=%s",
                    len(targets), len(target_tokens), target_type, job_id)

    meta = {
        "job_id": job_id,
        "status": "queued",
        "created_at": int(time.time()),
        "mirna_id": mirna_id,
        "mirna_seq": mirna_seq,
        "tools": tools,
        "workflow": workflow,
        "genome": genome,
        "cores": cores,
        "modifications": modifications,
        "shift": shift,
        "pre_id": pre_id,
        "targets": targets,
        "target_file": target_file,
    }
    _write_meta(job_id, meta)
    task = run_job.delay(job_id)
    meta["task_id"] = task.id
    _write_meta(job_id, meta)
    logger.info("job queued job_id=%s mirna_id=%s tools=%s genome=%s",
                job_id, mirna_id, tools, data.get("genome", "hg38"))
    return jsonify({"job_id": job_id, "task_id": task.id}), 202


@app.route("/api/v1/jobs/<job_id>", methods=["GET"])
def job_status(job_id):
    if not os.path.exists(_job_meta_path(job_id)):
        return jsonify({"error": "job not found"}), 404
    meta = _load_meta(job_id)
    progress = _load_progress(job_id)
    if progress is not None:
        meta["progress"] = progress
    resp = jsonify(meta)
    resp.headers["Cache-Control"] = "no-store"
    return resp


@app.route("/api/v1/jobs/<job_id>/result", methods=["GET"])
def job_result(job_id):
    if not os.path.exists(_job_meta_path(job_id)):
        return jsonify({"error": "job not found"}), 404
    meta = _load_meta(job_id)
    if meta.get("status") != "succeeded":
        logger.warning("result requested but job not succeeded job_id=%s status=%s",
                       job_id, meta.get("status"))
        return jsonify({"error": "job not completed", "status": meta.get("status")}), 409

    output_dir = meta.get("result_path")
    if not output_dir or not os.path.exists(output_dir):
        logger.error("result directory missing job_id=%s result_path=%s", job_id, output_dir)
        return jsonify({"error": "result not found"}), 404

    # Query parameters
    sort_by = request.args.get("sortBy", "tool_count")
    order   = request.args.get("ascendOrDescend", "desc")
    keyword = request.args.get("keyword", "").strip() or None
    try:
        offset = int(request.args.get("offset", 0))
        number = int(request.args.get("number", 20))
    except ValueError:
        return jsonify({"error": "offset and number must be integers"}), 400

    if sort_by not in ("gene_label", "tool_count"):
        return jsonify({"error": "sortBy must be gene_label or tool_count"}), 400
    if order.lower() not in ("asc", "desc", "ascend", "descend"):
        return jsonify({"error": "ascendOrDescend must be asc or desc"}), 400
    if offset < 0 or number < 1 or number > 1000:
        return jsonify({"error": "offset must be >= 0 and number must be between 1 and 1000"}), 400

    try:
        # lncRNA results are aggregated by transcript ID with no gene mapping
        # (see lncrna_results.py); the gene flow uses the RefSeq/symbol parser.
        if meta.get("workflow") == "mir-lncrna":
            db_path = ensure_lncrna_db(output_dir)
        else:
            db_path = ensure_db(output_dir)
        data = query_genes(db_path, sort_by=sort_by, order=order,
                           offset=offset, number=number, keyword=keyword)
    except Exception as e:
        logger.error("result processing failed job_id=%s error=%s", job_id, e)
        return jsonify({"error": "result processing failed: {}".format(str(e))}), 500

    data.update({
        "job_id":  job_id,
        "sort_by": sort_by,
        "order":   order,
        "offset":  offset,
        "number":  number,
        "keyword": keyword,
    })
    logger.info("result queried job_id=%s sort_by=%s order=%s offset=%d number=%d keyword=%s total=%d total_genes=%d",
                job_id, sort_by, order, offset, number, keyword, data["total"], data["total_genes"])
    return jsonify(data)


@app.route("/api/v1/jobs/<job_id>/network", methods=["GET"])
def job_network(job_id):
    """Tripartite gene <-> miRNA <-> lncRNA network for a mir-network job.

    Builds (and caches) network.db from both prediction pools, then shapes the
    graph: with user (gene, lncRNA) pairs it returns only the ceRNA co-target
    bridges; without pairs it returns a bounded discovery graph over the
    top-degree genes/lncRNAs (topGenes/topLncrna query params)."""
    if not os.path.exists(_job_meta_path(job_id)):
        return jsonify({"error": "job not found"}), 404
    meta = _load_meta(job_id)
    if meta.get("workflow") != NETWORK_WORKFLOW:
        return jsonify({"error": "network is only available for the mir-network workflow"}), 400
    if meta.get("status") != "succeeded":
        logger.warning("network requested but job not succeeded job_id=%s status=%s",
                       job_id, meta.get("status"))
        return jsonify({"error": "job not completed", "status": meta.get("status")}), 409

    output_dir = meta.get("result_path")
    if not output_dir or not os.path.exists(output_dir):
        logger.error("result directory missing job_id=%s result_path=%s", job_id, output_dir)
        return jsonify({"error": "result not found"}), 404

    try:
        top_genes = int(request.args.get("topGenes", DEFAULT_TOP_GENES))
        top_lncrna = int(request.args.get("topLncrna", DEFAULT_TOP_LNCRNA))
    except ValueError:
        return jsonify({"error": "topGenes and topLncrna must be integers"}), 400
    if top_genes < 1 or top_lncrna < 1:
        return jsonify({"error": "topGenes and topLncrna must be >= 1"}), 400

    pairs = load_pairs(_job_path(job_id))
    try:
        data = build_network(output_dir, pairs=pairs,
                             top_genes=top_genes, top_lncrna=top_lncrna)
    except Exception as e:
        logger.error("network build failed job_id=%s error=%s", job_id, e)
        return jsonify({"error": "network build failed: {}".format(str(e))}), 500

    data["job_id"] = job_id
    summary = data.get("summary", {})
    logger.info("network queried job_id=%s mode=%s genes=%s mirnas=%s lncrnas=%s edges=%s pairs=%s",
                job_id, data.get("mode"), summary.get("gene_count"), summary.get("mirna_count"),
                summary.get("lncrna_count"), summary.get("edge_count"), summary.get("pair_count"))
    resp = jsonify(data)
    resp.headers["Cache-Control"] = "no-store"
    return resp


@app.route("/api/v1/jobs/<job_id>/result/download", methods=["GET"])
def job_result_download(job_id):
    if not os.path.exists(_job_meta_path(job_id)):
        return jsonify({"error": "job not found"}), 404
    meta = _load_meta(job_id)
    if meta.get("status") != "succeeded":
        return jsonify({"error": "job not completed", "status": meta.get("status")}), 409

    result_dir = meta.get("result_path")
    if not result_dir or not os.path.exists(result_dir):
        return jsonify({"error": "result not found"}), 404

    archive_path = os.path.join(_job_path(job_id), "result.zip")
    if not os.path.exists(archive_path):
        # Build the archive ourselves (instead of shutil.make_archive) so we
        # can also include files that live in the job dir (sibling of
        # result_dir), not inside result_dir itself: targets.txt plus the run
        # inputs job.json and mirna.fa. Inside the zip, mirror result_dir/* at
        # the top level (matching shutil.make_archive's layout) and place those
        # sibling files alongside it.
        with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for root, _dirs, files in os.walk(result_dir):
                for name in files:
                    abs_path = os.path.join(root, name)
                    zf.write(abs_path, os.path.relpath(abs_path, result_dir))
            for sibling in ("targets.txt", "job.json", "mirna.fa"):
                sibling_path = os.path.join(_job_path(job_id), sibling)
                if os.path.exists(sibling_path):
                    zf.write(sibling_path, sibling)

    logger.info("result downloaded job_id=%s", job_id)
    return send_file(archive_path, as_attachment=True, download_name="{}_result.zip".format(job_id))


@app.route("/api/v1/jobs/<job_id>/enrichment", methods=["POST"])
def job_enrichment(job_id):
    if not os.path.exists(_job_meta_path(job_id)):
        return jsonify({"error": "job not found"}), 404

    data = request.get_json(force=True, silent=True) or {}
    genes = data.get("genes", [])
    organism = data.get("organism", "Human")
    cutoff = data.get("cutoff", 0.05)

    if not genes or not isinstance(genes, list):
        return jsonify({"error": "genes must be a non-empty list of gene symbols"}), 400

    try:
        cutoff = float(cutoff)
    except (TypeError, ValueError):
        return jsonify({"error": "cutoff must be a number"}), 400

    job_dir = _job_path(job_id)
    gene_list_csv = os.path.join(job_dir, "gene_list.csv")

    with open(gene_list_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_label"])
        for gene in genes:
            writer.writerow([gene])

    logger.info("enrichment started job_id=%s genes=%d organism=%s cutoff=%s",
                job_id, len(genes), organism, cutoff)

    # -u: unbuffered stdout/stderr so output is flushed as it happens.
    # stdout inherits gunicorn's stdout (visible in supervisord logs in real
    # time). Only stderr is captured so we can surface failure details — and
    # we cap it at 64KB so a runaway traceback can't fill memory.
    cmd = [
        "python3.6", "-u", "/opt/v2/enrichment_analysis.py",
        "-f", gene_list_csv,
        "-o", organism,
        "-c", str(cutoff),
        "-d", job_dir,
    ]
    try:
        subprocess.run(
            cmd,
            check=True,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as e:
        stderr_tail = (e.stderr or "")[-65536:].strip()
        logger.error(
            "enrichment failed job_id=%s rc=%s\nstderr=%s",
            job_id, e.returncode, stderr_tail,
        )
        return jsonify({
            "error": "enrichment analysis failed",
            "stderr": stderr_tail,
        }), 500

    logger.info("enrichment succeeded job_id=%s outdir=%s", job_id, job_dir)
    return jsonify({"job_id": job_id, "outdir": job_dir})


def _parse_enrichr_report(path):
    rows = []
    with open(path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            entry = {
                "term": row.get("Term"),
                "overlap": row.get("Overlap"),
                "p_value": _to_float(row.get("P-value")),
                "adjusted_p_value": _to_float(row.get("Adjusted P-value")),
                "odds_ratio": _to_float(row.get("Odds Ratio")),
                "combined_score": _to_float(row.get("Combined Score")),
                "genes": [g for g in (row.get("Genes") or "").split(";") if g],
            }
            rows.append(entry)
    return rows


def _to_float(value):
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


@app.route("/api/v1/jobs/<job_id>/enrichment", methods=["GET"])
def get_job_enrichment(job_id):
    if not os.path.exists(_job_meta_path(job_id)):
        return jsonify({"error": "job not found"}), 404

    job_dir = _job_path(job_id)
    gene_list_csv = os.path.join(job_dir, "gene_list.csv")
    if not os.path.exists(gene_list_csv):
        return jsonify({"error": "enrichment has not been run for this job"}), 404

    report_paths = sorted(glob.glob(os.path.join(job_dir, "*.enrichr.reports.txt")))

    databases = {}
    for path in report_paths:
        fname = os.path.basename(path)
        gene_set = fname.split(".enrichr.reports.txt")[0]
        if "." in gene_set:
            gene_set = gene_set.rsplit(".", 1)[0]
        try:
            databases[gene_set] = _parse_enrichr_report(path)
        except Exception as e:
            logger.error("failed to parse enrichment report job_id=%s file=%s error=%s",
                         job_id, fname, e)
            databases[gene_set] = []

    has_dotplot = os.path.exists(os.path.join(job_dir, "enrichment_dotplot.png"))
    logger.info("enrichment results queried job_id=%s databases=%d",
                job_id, len(databases))
    return jsonify({
        "job_id": job_id,
        "databases": databases,
        "has_dotplot": has_dotplot,
    })


@app.route("/api/v1/jobs/<job_id>/enrichment/dotplot", methods=["GET"])
def get_job_enrichment_dotplot(job_id):
    if not os.path.exists(_job_meta_path(job_id)):
        return jsonify({"error": "job not found"}), 404
    dotplot_path = os.path.join(_job_path(job_id), "enrichment_dotplot.png")
    if not os.path.exists(dotplot_path):
        return jsonify({"error": "dotplot not found"}), 404
    return send_file(dotplot_path, mimetype="image/png")


@app.route("/api/v1/jobs/<job_id>/kill", methods=["POST"])
def kill_job(job_id):
    if not os.path.exists(_job_meta_path(job_id)):
        return jsonify({"error": "job not found"}), 404
    meta = _load_meta(job_id)
    if meta.get("status") not in ("queued", "running"):
        logger.warning("kill requested but job not killable job_id=%s status=%s",
                       job_id, meta.get("status"))
        return jsonify({"error": "job is not killable", "status": meta.get("status")}), 409
    task_id = meta.get("task_id")
    if task_id:
        celery_app.control.revoke(task_id, terminate=True, signal="SIGTERM")
    meta["status"] = "killed"
    meta["finished_at"] = int(time.time())
    _write_meta(job_id, meta)
    logger.info("job killed job_id=%s", job_id)
    return jsonify({"job_id": job_id, "status": "killed"}), 200


@app.route("/api/v1/jobs/<job_id>", methods=["DELETE"])
def delete_job(job_id):
    if not os.path.exists(_job_meta_path(job_id)):
        return jsonify({"error": "job not found"}), 404
    shutil.rmtree(_job_path(job_id), ignore_errors=True)
    logger.info("job deleted job_id=%s", job_id)
    return jsonify({"status": "deleted"}), 200


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080)
