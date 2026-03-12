import json
import os
import shutil
import subprocess
import time
import uuid

from flask import Flask, jsonify, request, send_file

from app_v1.celery_app import celery_app
from app_v1.logger import get_logger

logger = get_logger()

BASE_DIR = os.environ.get("ISOTAR_JOB_DIR", "/opt/out/jobs")
os.makedirs(BASE_DIR, exist_ok=True)

app = Flask(__name__)

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

    if not tools or not mirna_id:
        logger.warning("job rejected: missing tools or mirna_id")
        return jsonify({"error": "tools and mirna_id are required"}), 400

    if not isinstance(modifications, list):
        logger.warning("job rejected: modifications not a list mirna_id=%s", mirna_id)
        return jsonify({"error": "modifications must be a list of strings"}), 400

    if shift is not None and not isinstance(shift, str):
        logger.warning("job rejected: invalid shift mirna_id=%s shift=%s", mirna_id, shift)
        return jsonify({"error": "shift must be a string in format 'left|right' (e.g. '-4|-6')"}), 400

    job_id = str(uuid.uuid4())
    job_dir = _job_path(job_id)
    os.makedirs(job_dir, exist_ok=True)

    meta = {
        "job_id": job_id,
        "status": "queued",
        "created_at": int(time.time()),
        "mirna_id": mirna_id,
        "tools": tools,
        "genome": data.get("genome", "hg38"),
        "cores": data.get("cores", 1),
        "modifications": modifications,
        "shift": shift,
        "pre_id": pre_id,
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

    result_dir = meta.get("result_path")
    if not result_dir or not os.path.exists(result_dir):
        logger.error("result directory missing job_id=%s result_path=%s", job_id, result_dir)
        return jsonify({"error": "result not found"}), 404

    archive_path = os.path.join(_job_path(job_id), "result.zip")
    if not os.path.exists(archive_path):
        shutil.make_archive(archive_path.replace(".zip", ""), "zip", result_dir)

    logger.info("result downloaded job_id=%s", job_id)
    return send_file(archive_path, as_attachment=True, download_name="{}_result.zip".format(job_id))


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
