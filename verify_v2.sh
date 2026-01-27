#!/usr/bin/env bash
set -euo pipefail

python /opt/v2/mirna_processing.py hsa-miR-495-3p -o /work/mirna_test.fa

python /opt/v2/mirna_predicting.py \
  -c 4 \
  -i /work/mirna_test.fa \
  -t RNAhybrid \
  -g hg38 \
  -o /work/mirna_predict_test

ls -la /work/mirna_test.fa
ls -la /work/mirna_predict_test/miRanda

echo "OK: functional test completed"
