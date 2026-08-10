# Changelog

All notable changes to the isoTar-v2 backend are documented here. This project
adheres to [Semantic Versioning](https://semver.org). Each version corresponds
to a `frankfeng78/isotar-v2` Docker image tag and a `vX.Y.Z` git tag.

New releases are cut with `scripts/release.sh <major|minor|patch>`.

## 0.3.22 - 2026-08-10

- fix(pita): give each PITA run a private cwd to stop cross-worker corruption

## 0.3.21 - 2026-08-04

- (no changes recorded since v0.3.20)

## 0.3.20 - 2026-08-04

- feat(download): include enrichment results in the job result zip

## 0.3.19 - 2026-08-03

- chore: commit stray trailing newlines in README
- fix(release): don't block a release on untracked files
- fix(mirmap): restore v1 parity — drop dG gate, fix seed lengths and multi-site parsing
- fix(predict): restore v1 parity for miRanda and PITA
- chore(release): v0.3.18
- feat(network): optionally restrict prediction to the ceRNA pair targets
- fix(result): read the requested pool for network jobs instead of 500ing
- fix(progress): accumulate per-tool elapsed instead of a start/finish bracket
- chore(release): v0.3.17
- feat(network): require a numeric score on each ceRNA (gene, lncRNA) pair
- chore(release): v0.3.16
- chore(release): v0.3.15
- fix(targetscan): sanitize filenames so variant miRNAs don't crash the tool
- feat(mir-network): support multiple sequence variants per miRNA
- feat(mir-network): support per-miRNA shift/modification variants as graph nodes
- chore(release): v0.3.14
- chore(release): v0.3.13
- feat(mir-lncrna): filter predictions by user-provided lncRNA targets
- test(targets): assert Ensembl IDs are rejected on the gene path
- feat(targets): add lncRNA/gene target validation endpoint
- chore(release): v0.3.12
- fix(network): serve merged per-pool progress for mir-network jobs
- chore(release): v0.3.11
- feat(mir-network): support per-miRNA precursor selection via pre_ids
- chore(release): v0.3.10
- fix the error that -s parsed incorrectly for mirna_processing.py
- chore(release): v0.3.9
- chore(release): v0.3.8
- chore(release): v0.3.7
- feat(mirna): support custom user-supplied miRNA sequences
- chore(mirmap): remove miRmap 1.x from the image
- fix(mirmap): stop PYTHONPATH=/opt/miRmap/src shadowing miRmap 2
- feat(mirmap): migrate miRmap to miRmap 2 on python3.11
- chore(release): v0.3.6
- feat(network): add mir-network workflow for gene↔miRNA↔lncRNA visualization
- fix header-row fix for PITA
- chore(release): v0.3.5
- fix(lncrna): make miRmap resilient to oversized transcripts and partial failures
- chore(release): v0.3.4
- feat(lncrna): support PITA in the miR-LncRNA workflow
- chore(release): v0.3.3
- feat(download): include job.json and mirna.fa in result.zip
- chore(release): v0.3.2
- fix(lncrna): extract non-human Ensembl/FlyBase/WormBase transcript IDs
- chore(release): v0.3.1
- feat(lncrna): add miR-LncRNA target prediction pipeline
- chore(release): v0.3.0
- chore: harden Dockerfile for species-aware metadata + gitignore local opt mirrors
- make mirna_processing.py support other species
- add mirna meta json files for other species
- Replace non-ASCII arrow in logger.py docstring with ASCII ->
- Make RNAhybrid distribution set and KEGG library species-aware
- introduce versioning (VERSION file, release script, /api/v1/version)
- fix enrichment 500: self-contained logging, drop broken logger import
- fix per-tool progress refresh and 4+-tool overlap views
- fix prediction parsers (miRanda/miRmap) and include targets.txt in result zip
- Integrate DMISO perfect-seed-match filter into parse_result
- fix docker logs error
- minor changes
- fix an encoding bug and add pre-commit hooks
- minor changes
- increase timeout to avoid enrichment fail
- add per-job and global core limits
- add target-based filtering and parallel TargetScan support
- fix the Targetscan result file name error
- control concurrency job number
- support target file
- implement enrichment feature
- update endpoint GET /api/v1/jobs/<job_id>/result
- enable enrichment_analysis to specify output directory
- enable enrichment_analysis to support to customize organism, cutoff and gene list
- add README.md
- support to map gene id to gene label when getting result
- move reference_mapping.db to app_v1 folder
- add reference mapping files
- add more species reference files
- fix a keyerror bug
- fix a bug
- add endpoint /api/v1/<job_id>/result
- add logging and remove sanitizing in mirna_predicting.py
- make parses_result to support Ensembl transcript ID (like: ENST00000284637.14) and RefSeq mRNA ID (like: NM_001164664)
- add scripts
- support output progress
- update endpoints
- add missed script file
- fix the issue that app_v1 server is not run
- add AGENTS.md and CLAUDE.md
- minor bugs fix
- add flask app v1
- make single Dockerfile work
- add tools but not include 3utr in TargetScan and CLAPACK due to size
- update Dockerfile to integrate 6 tools into single image
- update core Dockerfile and dmiso Dockerfile
- divide tools to three separate docker images
- divide dmiso and mirmap to separate Docker file
- add h5py==2.10.0 for dmiso
- fix path issue
- add v2 python files
- integrate DMISO
- add Dockerfile
- Initial commit

## 0.3.18 - 2026-07-16

- feat(network): optionally restrict prediction to the ceRNA pair targets
- fix(result): read the requested pool for network jobs instead of 500ing
- fix(progress): accumulate per-tool elapsed instead of a start/finish bracket

## 0.3.17 - 2026-07-15

- feat(network): require a numeric score on each ceRNA (gene, lncRNA) pair

## 0.3.16 - 2026-07-06

- (no changes recorded since v0.3.15)

## 0.3.15 - 2026-07-06

- fix(targetscan): sanitize filenames so variant miRNAs don't crash the tool
- feat(mir-network): support multiple sequence variants per miRNA
- feat(mir-network): support per-miRNA shift/modification variants as graph nodes

## 0.3.14 - 2026-07-04

- (no changes recorded since v0.3.13)

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
