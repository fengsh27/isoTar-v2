#!/usr/bin/env bash
set -euo pipefail

CORE_IMAGE=${CORE_IMAGE:-isotar-core:16.04}
MIRMAP_IMAGE=${MIRMAP_IMAGE:-isotar-mirmap:3.11}
DMISO_IMAGE=${DMISO_IMAGE:-isotar-dmiso:3.6}
WORK_DIR=${WORK_DIR:-"$(pwd)/work"}

mkdir -p "${WORK_DIR}"

echo "Generating miRNA FASTA in core image..."
docker run --rm \
  -v "${WORK_DIR}":/work \
  "${CORE_IMAGE}" \
  python /opt/v2/mirna_processing.py hsa-miR-495-3p -o /work/mirna_test.fa

echo "Running miRanda (core image)..."
docker run --rm \
  -v "${WORK_DIR}":/work \
  "${CORE_IMAGE}" \
  python /opt/v2/mirna_predicting.py \
    -c 1 \
    -i /work/mirna_test.fa \
    -t miRanda \
    -g hg38 \
    -o /work/mirna_predict_miranda

echo "Running miRmap (mirmap image)..."
docker run --rm \
  -v "${WORK_DIR}":/work \
  "${MIRMAP_IMAGE}" \
  python /opt/v2/mirna_predicting.py \
    -c 1 \
    -i /work/mirna_test.fa \
    -t miRmap \
    -g hg38 \
    -o /work/mirna_predict_mirmap

echo "Running DMISO (dmiso image)..."
docker run --rm \
  -v "${WORK_DIR}":/work \
  "${DMISO_IMAGE}" \
  python /opt/v2/mirna_predicting.py \
    -c 1 \
    -i /work/mirna_test.fa \
    -t DMISO \
    -g hg38 \
    -o /work/mirna_predict_dmiso

echo "Checking outputs..."
ls -la /work/mirna_test.fa
ls -la /work/mirna_predict_miranda/miRanda || true
ls -la /work/mirna_predict_mirmap/miRmap || true
ls -la /work/mirna_predict_dmiso/DMISO || true

echo "OK: tools verification completed"
