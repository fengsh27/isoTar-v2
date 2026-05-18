import csv
import glob
import json
import os
import shutil
import subprocess
import time
import uuid

from flask import Flask, jsonify, request, send_file
from flask_cors import CORS

from app_v1.celery_app import celery_app
from app_v1.limits import MAX_CORES_PER_JOB, validate_cores
from app_v1.logger import get_logger
from app_v1.result_db import ensure_db, query_genes
from app_v1.target_resolver import (
    genome_to_species as _genome_to_species,
    resolve_targets as _resolve_targets,
)

logger = get_logger()

BASE_DIR = os.environ.get("ISOTAR_JOB_DIR", "/opt/out/jobs")
os.makedirs(BASE_DIR, exist_ok=True)

app = Flask(__name__)

_cors_origins = os.environ.get("ISOTAR_CORS_ORIGINS", "*")
if _cors_origins != "*":
    _cors_origins = [o.strip() for o in _cors_origins.split(",") if o.strip()]
CORS(app, resources={r"/api/*": {"origins": _cors_origins}})


def _job_path(job_id):
    return os.path.join(BASE_DIR, job_id)


def _job_meta_path(job_id):
    return os.path.join(_job_path(job_id), "job.json")


def _progress_path(job_id):
    return os.path.join(_job_path(job_id), "output", "progress.json")


def _write_meta(job_id, data):
    with open(_job_meta_path(job_id), "w") as f:
        json.dump(data, f, indent=2)


def _load_meta(job_id):
    with open(_job_meta_path(job_id), "r") as f:
        return json.load(f)


def _load_progress(job_id):
    path = _progress_path(job_id)
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r") as f:
            return json.load(f)
    except (ValueError, IOError):
        return None


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
        mirna_id = meta["mirna_id"]
        tools = meta["tools"]
        genome = meta.get("genome", "hg38")
        cores = int(meta.get("cores", 1))
        modifications = meta.get("modifications", [])
        shift = meta.get("shift")
        pre_id = meta.get("pre_id")
        target_file = meta.get("target_file")

        logger.info("job started job_id=%s mirna_id=%s tools=%s genome=%s cores=%s",
                    job_id, mirna_id, tools, genome, cores)

        fasta_path = os.path.join(job_dir, "mirna.fa")
        cmd = [
            "python2.7",
            "/opt/v2/mirna_processing.py",
            mirna_id,
            "-o",
            fasta_path,
        ]
        for mod in modifications:
            cmd.extend(["-m", mod])
        if shift:
            cmd.extend(["-s", shift])
        if modifications and shift:
            cmd.append("-b")
        if pre_id:
            cmd.extend(["--pre-id", pre_id])
        logger.info("step=processing job_id=%s", job_id)
        subprocess.check_call(cmd)
        meta["step"] = "processing"
        _write_meta(job_id, meta)

        meta["step"] = "predicting"
        _write_meta(job_id, meta)
        logger.info("step=predicting job_id=%s", job_id)

        other_tools = [t for t in tools if t != "miRmap"]
        if other_tools:
            logger.info("running python3.6 tools=%s job_id=%s", other_tools, job_id)
            cmd = [
                "python3.6",
                "/opt/v2/mirna_predicting.py",
                "-c", str(cores),
                "-i", fasta_path,
                "-t",
            ] + other_tools + ["-g", genome, "-o", output_dir]
            if target_file:
                cmd += ["-tf", target_file]
            subprocess.check_call(cmd)

        if "miRmap" in tools:
            logger.info("running python2.7 tools=['miRmap'] job_id=%s", job_id)
            cmd = [
                "python2.7",
                "/opt/v2/mirna_predicting.py",
                "-c", str(cores),
                "-i", fasta_path,
                "-t", "miRmap",
                "-g", genome,
                "-o", output_dir,
            ]
            if target_file:
                cmd += ["-tf", target_file]
            subprocess.check_call(cmd)

        meta["status"] = "succeeded"
        meta["finished_at"] = int(time.time())
        meta["result_path"] = output_dir
        _write_meta(job_id, meta)
        logger.info("job succeeded job_id=%s", job_id)
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


@app.route("/api/v1/jobs", methods=["POST"])
def submit_job():
    data = request.get_json(force=True, silent=True) or {}
    tools = data.get("tools", [])
    mirna_id = data.get("mirna_id")
    modifications = data.get("modifications", [])
    shift = data.get("shift")
    pre_id = data.get("pre_id")
    targets = data.get("targets")

    if not tools or not mirna_id:
        logger.warning("job rejected: missing tools or mirna_id")
        return jsonify({"error": "tools and mirna_id are required"}), 400

    if not isinstance(modifications, list):
        logger.warning("job rejected: modifications not a list mirna_id=%s", mirna_id)
        return jsonify({"error": "modifications must be a list of strings"}), 400

    if shift is not None and not isinstance(shift, str):
        logger.warning("job rejected: invalid shift mirna_id=%s shift=%s", mirna_id, shift)
        return jsonify({"error": "shift must be a string in format 'left|right' (e.g. '-4|-6')"}), 400

    if targets is not None:
        if not isinstance(targets, list) or not targets or not all(isinstance(g, str) for g in targets):
            return jsonify({"error": "targets must be a non-empty list of gene symbol strings"}), 400

    cores, cores_err = validate_cores(data.get("cores"))
    if cores_err:
        msg, status = cores_err
        logger.warning("job rejected: %s mirna_id=%s cores=%r", msg, mirna_id, data.get("cores"))
        return jsonify({"error": msg, "max_cores_per_job": MAX_CORES_PER_JOB}), status

    genome = data.get("genome", "hg38")

    # Resolve targets before creating any job state so a 400 leaves no orphan dir.
    accessions = None
    if targets:
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

    job_id = str(uuid.uuid4())
    job_dir = _job_path(job_id)
    os.makedirs(job_dir, exist_ok=True)

    target_file = None
    if accessions:
        target_file = os.path.join(job_dir, "targets.txt")
        with open(target_file, "w") as f:
            for acc in sorted(accessions):
                f.write(acc + "\n")
        logger.info("resolved %d target tokens to %d UTR identifiers job_id=%s",
                    len(targets), len(accessions), job_id)

    meta = {
        "job_id": job_id,
        "status": "queued",
        "created_at": int(time.time()),
        "mirna_id": mirna_id,
        "tools": tools,
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
        shutil.make_archive(archive_path.replace(".zip", ""), "zip", result_dir)

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
    cmd = [
        "python3.6", "-u", "/opt/v2/enrichment_analysis.py",
        "-f", gene_list_csv,
        "-o", organism,
        "-c", str(cutoff),
        "-d", job_dir,
    ]
    try:
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as e:
        logger.error(
            "enrichment failed job_id=%s rc=%s\nstdout=%s\nstderr=%s",
            job_id, e.returncode, e.stdout, e.stderr,
        )
        return jsonify({
            "error": "enrichment analysis failed",
            "stderr": (e.stderr or "").strip(),
        }), 500

    if result.stdout:
        logger.info("enrichment stdout job_id=%s\n%s", job_id, result.stdout)
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
