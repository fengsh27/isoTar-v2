# isoTar-v2

A microRNA (miRNA) target prediction platform that runs six computational tools in parallel and aggregates results into a searchable database. Supports canonical miRNA sequences as well as isomiR variants (modified, shifted, or truncated).

## Overview

isoTar-v2 accepts a miRNA identifier plus optional sequence modifications, generates the corresponding FASTA sequences, and dispatches prediction jobs to any combination of six tools. Jobs are processed asynchronously via Celery. Results are stored in SQLite and exposed through a REST API with pagination, sorting, and gene label filtering.

## Architecture

```
Client
  └── HTTP :8080
        └── Nginx (reverse proxy)
              ├── /api/v1/*  →  Gunicorn :5001  (app_v1 — Python 3.6)
              └── /*         →  Gunicorn :5000  (legacy app — Python 2.7)

app_v1 (Flask + Celery)
  └── run_job task
        ├── python2.7 mirna_processing.py   (FASTA generation)
        ├── python3.6 mirna_predicting.py   (all tools except miRmap)
        └── python2.7 mirna_predicting.py   (miRmap only)

Supervisor manages: nginx, rabbitmq, gunicorn (legacy), app_v1 gunicorn, celery worker
```

## Prediction Tools

| Tool | Version | Method |
|------|---------|--------|
| miRanda | 3.3 | Seed-based thermodynamic scoring |
| miRmap | 1.1 | Comprehensive thermodynamic model (Python 2.7) |
| RNAhybrid | 2.1.2 | RNA–RNA interaction energy |
| PITA | v6 | Position-and-context-dependent ddG scoring |
| TargetScan | 7.0 | Seed match + conservation scoring |
| DMISO | — | Deep learning (TensorFlow/Keras) |

## Supported Genomes / Species

| Code | Species |
|------|---------|
| `hg19` | Human (GRCh37) |
| `hg38` | Human (GRCh38) |
| `mmu` | House mouse (GRCm38) |
| `rno` | Norway rat (RGSC6/rn6) |
| `dre` | Zebrafish (GRCz11) |
| `dme` | Fruit fly (Release 6) |
| `cel` | Roundworm (WBcel235) |
| `cfa` | Dog (CanFam3.1) |
| `mdo` | Gray short-tailed opossum (MonDom5) |
| `mml` | Rhesus macaque (Mmul_8.0.1) |
| `ptr` | Chimpanzee (Pan_tro3.0) |

## Docker Build

```bash
# Step 1 — base image (Ubuntu 16.04, Python 3.6 compiled from source, build tools)
# Only needed once or when base dependencies change
docker build -t frankfeng78/isotar-v2-base:0.2.x -f isotar-base.Dockerfile .

# Step 2 — final production image
docker build -t frankfeng78/isotar-v2:0.2.x -f Dockerfile .
```

## Running

```bash
docker run -d \
  --name isotar \
  -p 8080:8080 \
  -v /path/to/jobs:/opt/out/jobs \
  frankfeng78/isotar-v2:0.2.x
```

The API is available at `http://localhost:8080/api/v1/`.

## REST API

### Submit a job

```
POST /api/v1/jobs
```

Request body:

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `mirna_id` | string | yes | miRNA identifier (e.g. `hsa-miR-21-5p`) |
| `tools` | string[] | yes | Tools to run: `miRanda`, `miRmap`, `RNAhybrid`, `PITA`, `Targetscan`, `DMISO` |
| `genome` | string | yes | Species code (see table above) |
| `cores` | int | no | Number of CPU cores (default: 1) |
| `modifications` | string[] | no | Nucleotide modifications |
| `shift` | string | no | Sequence shift (e.g. `"-2|-3"`) |
| `pre_id` | string | no | Precursor miRNA ID |

Response `202`:
```json
{ "job_id": "<uuid>", "task_id": "<celery-task-id>" }
```

---

### Get job status

```
GET /api/v1/jobs/<job_id>
```

Response includes `status` (`queued` | `running` | `succeeded` | `failed` | `killed`) and a `progress` object showing per-tool status when the prediction step is active.

---

### Query results

```
GET /api/v1/jobs/<job_id>/result
```

Query parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sortBy` | `tool_count` | Sort by `tool_count` or `gene_label` |
| `ascendOrDescend` | `desc` | `asc` or `desc` |
| `offset` | `0` | Pagination offset |
| `number` | `20` | Page size (max 1000) |
| `geneLabel` | — | Filter by gene symbol (substring match) |

Response:
```json
{
  "total_genes": 1500,
  "total": 42,
  "genes": [
    {
      "gene_id": "NM_001234",
      "gene_label": "TP53",
      "gene_name": "tumor protein p53",
      "tool_count": 4,
      "tools": ["DMISO", "PITA", "RNAhybrid", "Targetscan"]
    }
  ],
  "venn": {
    "sets": { "PITA": 310, "RNAhybrid": 220, "Targetscan": 180 },
    "intersections": { "PITA&RNAhybrid": 85, "PITA&Targetscan": 60 }
  }
}
```

---

### Download results

```
GET /api/v1/jobs/<job_id>/result/download
```

Returns a ZIP archive of all raw prediction output files.

---

### Kill a job

```
POST /api/v1/jobs/<job_id>/kill
```

Terminates a queued or running job.

---

### Delete a job

```
DELETE /api/v1/jobs/<job_id>
```

Removes all files and metadata for the job.

---

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `ISOTAR_JOB_DIR` | `/opt/out/jobs` | Root directory for job storage |
| `ISOTAR_LOG_DIR` | `/opt/out` | Log file directory |
| `ISOTAR_LOG_LEVEL` | `INFO` | Logging level |
| `ISOTAR_REFERENCE_MAPPING_DB` | `<app_v1>/reference_mapping.db` | Gene ID → symbol/name mapping database |
| `CELERY_BROKER_URL` | `amqp://` | RabbitMQ broker URL |
| `CELERY_RESULT_BACKEND` | `rpc://` | Celery result backend |

## Project Structure

```
isoTar-v2/
├── app/                        # Legacy Flask application (Python 2.7)
├── app_v1/                     # REST API v1 (Python 3.6)
│   ├── app.py                  # Flask app and Celery tasks
│   ├── celery_app.py           # Celery configuration
│   ├── result_db.py            # SQLite result database
│   ├── parse_result.py         # Tool output parsers
│   ├── logger.py               # Date-rotating file logger
│   └── reference_mapping.db   # Gene ID → label/name mapping
├── v2/                         # Prediction engine scripts
│   ├── mirna_processing.py     # miRNA FASTA generation
│   ├── mirna_predicting.py     # Tool orchestration and parallelism
│   ├── parse_result.py         # Result parsing library
│   ├── result_db.py            # Result database builder
│   └── enrichment_analysis.py # Gene enrichment analysis
├── v2/opt/reference_files/     # 3' UTR FASTA files per species
├── Dockerfile                  # Production image (builds on base)
├── isotar-base.Dockerfile      # Base image (Ubuntu 16.04 + Python 3.6)
├── app_supervisord.conf        # Supervisor process config
├── nginx.conf / app.conf       # Nginx configuration
└── CLAUDE.md                   # Developer notes
```
