import json
import os
import shutil
import subprocess
import time
import uuid

from flask import Flask, jsonify, request, send_file

from app_v1.celery_app import celery_app

BASE_DIR = os.environ.get("ISOTAR_JOB_DIR", "/opt/out/jobs")
os.makedirs(BASE_DIR, exist_ok=True)

app = Flask(__name__)

def _job_path(job_id):
    return os.path.join(BASE_DIR, job_id)


def _job_meta_path(job_id):
    return os.path.join(_job_path(job_id), "job.json")


def _write_meta(job_id, data):
    with open(_job_meta_path(job_id), "w") as f:
        json.dump(data, f, indent=2)


def _load_meta(job_id):
    with open(_job_meta_path(job_id), "r") as f:
        return json.load(f)


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
        operations = meta.get("operations", [])
        tools = meta["tools"]
        genome = meta.get("genome", "hg38")
        cores = int(meta.get("cores", 1))

        fasta_path = os.path.join(job_dir, "mirna.fa")
        cmd = [
            "python2.7",
            "/opt/v2/mirna_processing.py",
            mirna_id,
            "-o",
            fasta_path,
        ]
        for op in operations:
            cmd.extend(op)
        subprocess.check_call(cmd)

        cmd = [
            "python3.6",
            "/opt/v2/mirna_predicting.py",
            "-c",
            str(cores),
            "-i",
            fasta_path,
            "-t",
        ] + tools + [
            "-g",
            genome,
            "-o",
            output_dir,
        ]
        subprocess.check_call(cmd)

        meta["status"] = "succeeded"
        meta["finished_at"] = int(time.time())
        meta["result_path"] = output_dir
        _write_meta(job_id, meta)
    except subprocess.CalledProcessError as e:
        meta["status"] = "failed"
        meta["finished_at"] = int(time.time())
        meta["error"] = "Command failed: {}".format(e)
        _write_meta(job_id, meta)
    except Exception as e:
        meta["status"] = "failed"
        meta["finished_at"] = int(time.time())
        meta["error"] = str(e)
        _write_meta(job_id, meta)


@app.route("/api/v1/jobs", methods=["POST"])
def submit_job():
    data = request.get_json(force=True, silent=True) or {}
    tools = data.get("tools", [])
    mirna_id = data.get("mirna_id")
    operations = data.get("operations", [])

    if not tools or not mirna_id:
        return jsonify({"error": "tools and mirna_id are required"}), 400

    job_id = str(uuid.uuid4())
    job_dir = _job_path(job_id)
    os.makedirs(job_dir, exist_ok=True)

    meta = {
        "job_id": job_id,
        "status": "queued",
        "created_at": int(time.time()),
        "tools": tools,
        "mirna_id": mirna_id,
        "operations": operations,
        "genome": data.get("genome", "hg38"),
        "cores": data.get("cores", 1),
    }
    _write_meta(job_id, meta)
    task = run_job.delay(job_id)
    return jsonify({"job_id": job_id, "task_id": task.id}), 202


@app.route("/api/v1/jobs/<job_id>", methods=["GET"])
def job_status(job_id):
    if not os.path.exists(_job_meta_path(job_id)):
        return jsonify({"error": "job not found"}), 404
    meta = _load_meta(job_id)
    resp = jsonify(meta)
    resp.headers["Cache-Control"] = "no-store"
    return resp


@app.route("/api/v1/jobs/<job_id>/result", methods=["GET"])
def job_result(job_id):
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

    return send_file(archive_path, as_attachment=True, download_name="{}_result.zip".format(job_id))


@app.route("/api/v1/jobs/<job_id>", methods=["DELETE"])
def delete_job(job_id):
    if not os.path.exists(_job_meta_path(job_id)):
        return jsonify({"error": "job not found"}), 404
    shutil.rmtree(_job_path(job_id), ignore_errors=True)
    return jsonify({"status": "deleted"}), 200


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080)
