# Running isoTar on your own server

isoTar is published as two Docker images, so a lab can run the complete tool
on its own machine instead of the public site at https://isotar.osumc.edu:

| Image | What it is | Size |
|-------|------------|------|
| [`frankfeng78/isotar-v2`](https://hub.docker.com/r/frankfeng78/isotar-v2) | The prediction tools (miRanda, miRmap, RNAhybrid, PITA, TargetScan, DMISO), the job queue and the REST API | about 2.6 GB download |
| [`frankfeng78/isotar-v2-frontend`](https://hub.docker.com/r/frankfeng78/isotar-v2-frontend) | The web interface | about 110 MB download |

Reasons to run it locally:

- **More CPU per job.** The public server caps a job at 8 cores and runs four
  jobs at a time. On your own machine you set both limits, and the slow tools
  (RNAhybrid above all) scale with the cores you give them.
- **Your data stays with you.** Jobs, miRNA variants, target lists and results
  are written to a directory on your server and nowhere else.
- **No shared queue.** Your jobs never wait behind anyone else's.

There are two ways to use the images. **Option A** runs the website locally,
identical to the public one. **Option B** runs only the backend and drives it
from the command line or a script, which suits batch runs and pipelines. Both
produce the same results and the same result files.

## Requirements

- A Linux server (or macOS for small tests) with Docker and Docker Compose
  2.24 or newer. Run `docker compose version` to check.
- CPU: the tools are CPU-bound. 8 cores is enough to run one job at the public
  server's speed; the defaults below assume 32 cores can be given to isoTar.
- Disk: about 3 GB for the images plus the jobs you keep. A genome-wide job
  writes raw output for every tool; keep the data directory on a volume with
  tens of gigabytes free and delete jobs you no longer need.
- No GPU and no internet access are needed to run predictions. Internet is used
  only to pull the images and, if you enable them, for the optional PubMed and
  model-review features described below.

## Option A: the full website

1. Create a directory for the installation and download the two configuration
   files from this repository:

   ```bash
   mkdir isotar && cd isotar
   curl -LO https://raw.githubusercontent.com/fengsh27/isoTar-v2/main/docker-compose.yml
   curl -LO https://raw.githubusercontent.com/fengsh27/isoTar-v2/main/.env.template
   cp .env.template .env
   ```

   `.env` is optional. Leave it as it is for a first run; the settings it
   controls are explained under [Tuning](#tuning-for-your-machine) and
   [Optional features](#optional-features-and-credentials).

2. Start the stack:

   ```bash
   docker compose up -d
   ```

   The first start downloads the images. After that, starting takes a few
   seconds.

3. Open http://localhost:9501 (or the server's name instead of `localhost`).
   The interface and workflows are the same as the public site; the user guide
   under **Docs** applies unchanged.

4. Check the versions that are running:

   ```bash
   docker compose ps
   docker compose exec backend curl -s http://localhost:8080/api/v1/version
   ```

Everything isoTar writes goes to `./isotar-data` next to the compose file:
`jobs/<job_id>/` holds each job's inputs, raw tool outputs and result files.
Set `ISOTAR_DATA_DIR` in `.env` to use a different path, for example a large
data volume.

To stop isoTar, run `docker compose down`. Jobs are kept; `docker compose up -d`
brings them back.

## Option B: backend only

The backend image contains everything needed to run predictions and exposes
the same REST API the website uses, on port 8080.

```bash
docker run -d --name isotar \
  -p 8080:8080 \
  -e ISOTAR_DEFAULT_CORES=16 -e ISOTAR_MAX_CORES_PER_JOB=16 -e ISOTAR_MAX_CONCURRENT_JOBS=2 \
  -v "$PWD/isotar-data":/opt/out \
  frankfeng78/isotar-v2:0.3.35
```

A complete run from the shell:

```bash
API=http://localhost:8080/api/v1

# 1. Submit a job: hsa-miR-21-5p against human hg38 with four tools, 16 cores.
JOB=$(curl -s -X POST $API/jobs -H 'Content-Type: application/json' -d '{
  "workflow": "mir-target",
  "genome": "hg38",
  "mirna_id": "hsa-miR-21-5p",
  "tools": ["miRanda", "RNAhybrid", "PITA", "Targetscan"],
  "cores": 16
}' | python3 -c 'import sys, json; print(json.load(sys.stdin)["job_id"])')

# 2. Poll until status is "succeeded" (or "failed").
curl -s $API/jobs/$JOB | python3 -m json.tool

# 3. Browse the parsed results: targets sorted by how many tools agree.
curl -s "$API/jobs/$JOB/result?number=20" | python3 -m json.tool

# 4. Download everything as a zip: raw tool outputs, result.db,
#    predictions.csv, and the inputs that produced them.
curl -o isotar_$JOB.zip $API/jobs/$JOB/result/download
```

The same request body accepts a shift (`"shift": "-1|0"`), modifications
(`"modifications": ["8:A|G"]`), a target list (`"targets": ["PTEN", "PDCD4"]`),
the `mir-lncrna` workflow, and the `mir-network` workflow with its ceRNA pairs.
Every field, every endpoint and the response formats are documented in the
[REST API section of the README](README.md#rest-api).

## Tuning for your machine

Three settings decide how much of the machine isoTar uses. Set them in `.env`
(Option A) or with `-e` (Option B):

| Setting | Default | Meaning |
|---------|---------|---------|
| `ISOTAR_DEFAULT_CORES` | 8 | Cores a job gets when the submission does not say |
| `ISOTAR_MAX_CORES_PER_JOB` | 8 | Largest `cores` value a submission may ask for |
| `ISOTAR_MAX_CONCURRENT_JOBS` | 4 | Jobs running at the same time; further jobs queue |
| `ISOTAR_CPU_LIMIT` | 32 | Hard ceiling on the container (Option A only) |

Keep jobs times cores at or below the CPU limit. Two examples:

- **A 64-core server shared with other work.** `ISOTAR_DEFAULT_CORES=16`,
  `ISOTAR_MAX_CORES_PER_JOB=16`, `ISOTAR_MAX_CONCURRENT_JOBS=3`,
  `ISOTAR_CPU_LIMIT=48`.
- **A 16-core workstation for one user.** `ISOTAR_DEFAULT_CORES=12`,
  `ISOTAR_MAX_CORES_PER_JOB=12`, `ISOTAR_MAX_CONCURRENT_JOBS=1`,
  `ISOTAR_CPU_LIMIT=12`.

RNAhybrid takes most of a job's time, typically two thirds or more, and it is
the tool that benefits most from extra cores. If runtime matters more than
coverage, leave RNAhybrid out of the tool selection or, for network jobs,
submit with `restrict_to_pairs` so the tools scan only the pairs you supplied.

## Optional features and credentials

Predictions need no account or key of any kind. Two features of the
**miR-Network** workflow call external services and are switched on by
putting your own credentials in `.env`:

| Feature | Variables | Where to get them | When blank |
|---------|-----------|-------------------|------------|
| Literature verification (PubMed search for each predicted triplet) | `NCBI_API_KEY`, `NCBI_EMAIL` | Free, from your NCBI account settings under "API Key Management" | The "Verify against literature" button is not shown |
| Model review of the candidate papers | `ANTHROPIC_BASE_URL`, `ANTHROPIC_AUTH_TOKEN`, `ISOTAR_VERIFICATION_MODEL` | Your institution's Anthropic-compatible gateway (LiteLLM) and a key for it | The review option is disabled; PubMed verification still works |

Nothing else is affected by leaving these blank. The images contain no
credentials of ours; each installation uses its own, and the PubMed answers
and model verdicts are cached under the data directory so repeated
verifications cost nothing.

`.env.template` lists every other variable the server reads, with its default.

## Upgrading

New versions are announced in [CHANGELOG.md](CHANGELOG.md). To upgrade an
Option A installation, download the current `docker-compose.yml`, which pins a
matching pair of backend and frontend tags, then:

```bash
docker compose pull
docker compose up -d
```

Jobs and results are kept across upgrades. When a release changes the result
database layout, the changelog says so; existing jobs are rebuilt the first
time they are opened, which for a large job can take a minute. To rebuild them
all in advance:

```bash
docker compose exec -w / backend python3.6 -m app_v1.result_indexes
```

For Option B, stop and remove the old container and start the new tag with the
same `-v` data directory.

## Troubleshooting

**The web page loads but every submission fails.** The frontend cannot reach
the backend. Run `docker compose ps`; both services must be `running`. Look at
`docker compose logs backend` for the error. The backend takes about 30 seconds
to come up after a cold start because it waits for its message broker.

**`bind: address already in use` on start.** Port 9501 is taken on the host.
Set `ISOTAR_PORT` in `.env` to a free port.

**A job is rejected with HTTP 400 and "cores".** The submission asked for more
cores than `ISOTAR_MAX_CORES_PER_JOB` allows. Raise the limit or lower the
request.

**Jobs stay `queued`.** `ISOTAR_MAX_CONCURRENT_JOBS` jobs are already running.
Wait, or raise the limit if the machine has spare cores.

**TargetScan is refused for my species.** TargetScan datasets exist for human,
mouse, fly, zebrafish and worm only. Submit without TargetScan for the other
species; the API answers 400 with that reason.

**Permission denied writing to the data directory.** The container runs as root
and creates `isotar-data` itself on first start. If you pre-created it as
another user, make it writable or point `ISOTAR_DATA_DIR` at a fresh path.

**The "Verify against literature" button is missing.** `NCBI_API_KEY` is not
set in `.env`, or `.env` is not next to the compose file. After editing `.env`,
run `docker compose up -d` again so the backend picks it up.

## Where things are stored

```
isotar-data/
  jobs/<job_id>/
    job.json            parameters, status, timings
    mirna.fa            the sequences that were analysed (wild type and variants)
    targets.txt         the target list, if one was given
    output/             raw output of every tool, result.db, predictions.csv
  verification/         cached PubMed answers and model verdicts (shared by all jobs)
  *.log                 server and worker logs
```

A job is a plain directory. Copying `jobs/<job_id>/` to another installation
of the same version is enough to open it there.
