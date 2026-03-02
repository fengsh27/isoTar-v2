# Agent Guide

## Scope
- Repo combines legacy Flask (Python 2.7) with tool runners in `/opt/v2`.
- Multiple Python versions are used (2.7 for legacy stack + miRmap 1.x, 3.6 for DMISO and some runners).

## Python Version Guidance
- `v2/*.py` scripts are made compatible with Python 2.7/3.5 (no f-strings, no `exist_ok=True`).
- miRmap 1.x is Python 2.7 oriented; use `python2.7` when invoking miRmap in `v2/mirna_predicting.py`.
- DMISO should be run with `python3.6 /opt/DMISO/DMISO-main/dmiso.py`.

## Containers & Docker
- Base image file: `isotar-base.Dockerfile` (slow apt installs, Python 3.6 build).
- Final image file: `Dockerfile` (tool setup, app setup).
- Build sequence:
  - `docker build -t frankfeng78/isotar-v2-base:0.2.x -f isotar-base.Dockerfile .`
  - `docker build -t frankfeng78/isotar-v2:0.2.x -f Dockerfile .`

## Paths & Data
- UTR datasets: `/opt/human/hg19/3utr.fa` and `/opt/human/hg38/3utr.fasta`
- miRNA metadata: `/opt/resources/mature_pre_mirna_ext.json`
- TargetScan datasets are mounted or copied to `/opt/TargetScan/Datasets`

## Job Execution
- `v2/mirna_processing.py` generates miRNA FASTA.
- `v2/mirna_predicting.py` runs tools based on `-t` selection.
- Parallel branch uses multiprocessing `pool.map` with wrapper functions.

## Common Issues
- DMISO wrapper must preserve args; call script directly if needed.
- miRmap 1.x requires `dendropy` and Python 2.7; Python 3 will fail on `string.maketrans`.
- `subprocess.run` isn’t available in Python 2; ensure shim exists where used.

## Expectations
- Prefer minimal changes, keep legacy behavior.
- Avoid upgrading libraries unless required for compatibility.
