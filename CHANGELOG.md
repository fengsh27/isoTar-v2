# Changelog

All notable changes to the isoTar-v2 backend are documented here. This project
adheres to [Semantic Versioning](https://semver.org). Each version corresponds
to a `frankfeng78/isotar-v2` Docker image tag and a `vX.Y.Z` git tag.

New releases are cut with `scripts/release.sh <major|minor|patch>`.

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
