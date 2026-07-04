# Changelog

All notable changes to the isoTar-v2 backend are documented here. This project
adheres to [Semantic Versioning](https://semver.org). Each version corresponds
to a `frankfeng78/isotar-v2` Docker image tag and a `vX.Y.Z` git tag.

New releases are cut with `scripts/release.sh <major|minor|patch>`.

## 0.3.13 - 2026-07-04

- feat(mir-lncrna): filter predictions by user-provided lncRNA targets
- test(targets): assert Ensembl IDs are rejected on the gene path
- feat(targets): add lncRNA/gene target validation endpoint

## 0.3.12 - 2026-06-29

- fix(network): serve merged per-pool progress for mir-network jobs

## 0.3.11 - 2026-06-29

- feat(mir-network): support per-miRNA precursor selection via pre_ids

## 0.3.10 - 2026-06-26

- fix the error that -s parsed incorrectly for mirna_processing.py

## 0.3.9 - 2026-06-25

- (no changes recorded since v0.3.8)

## 0.3.8 - 2026-06-25

- (no changes recorded since v0.3.7)

## 0.3.7 - 2026-06-25

- feat(mirna): support custom user-supplied miRNA sequences
- chore(mirmap): remove miRmap 1.x from the image
- fix(mirmap): stop PYTHONPATH=/opt/miRmap/src shadowing miRmap 2
- feat(mirmap): migrate miRmap to miRmap 2 on python3.11

## 0.3.6 - 2026-06-23

- feat(network): add mir-network workflow for gene↔miRNA↔lncRNA visualization
- fix header-row fix for PITA

## 0.3.5 - 2026-06-22

- fix(lncrna): make miRmap resilient to oversized transcripts and partial failures

## 0.3.4 - 2026-06-22

- feat(lncrna): support PITA in the miR-LncRNA workflow

## 0.3.3 - 2026-06-16

- feat(download): include job.json and mirna.fa in result.zip

## 0.3.2 - 2026-06-16

- fix(lncrna): extract non-human Ensembl/FlyBase/WormBase transcript IDs

## 0.3.1 - 2026-06-16

- feat(lncrna): add miR-LncRNA target prediction pipeline

## 0.3.0 - 2026-06-09

- chore: harden Dockerfile for species-aware metadata + gitignore local opt mirrors
- make mirna_processing.py support other species
- add mirna meta json files for other species
- Replace non-ASCII arrow in logger.py docstring with ASCII ->
- Make RNAhybrid distribution set and KEGG library species-aware

## 0.2.27 - 2026-05-25

Baseline release — versioning introduced (`VERSION` file, git tags, this
changelog, `/api/v1/version` endpoint, version-stamped image). Notable changes
folded into this baseline:

- **Enrichment:** fix HTTP 500 caused by a broken `from logger import get_logger`
  (the module is not on the subprocess import path); use self-contained logging.
- **Result overlap:** `_venn_stats` now emits exclusive `combinations` and
  per-degree `degrees`, enabling the 4+-tool UpSet plot and consensus histogram.
- **Per-tool progress:** `_write_progress` merges with the existing
  `progress.json` so the second prediction phase no longer overwrites the first
  (all tools persist after a page refresh).
- **Prediction parsers:** fix miRanda parsing (handle `-quiet` summary lines)
  and miRmap header extraction (use the full RefSeq/Ensembl accession).
- **Result download:** include `targets.txt` in the result `.zip`.
- **Concurrency:** per-job and global core limits (`app_v1/limits.py`).
- **Tooling:** source-encoding pre-commit guard; regression tests.
