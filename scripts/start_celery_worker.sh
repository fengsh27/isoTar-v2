#!/bin/sh
# Launch the app_v1 celery worker with concurrency derived from MAX_CORE_PER_JOB.
# Concurrency = max(1, host_cores // MAX_CORE_PER_JOB).
# MAX_CORE_PER_JOB defaults to host_cores (i.e. concurrency = 1).

set -e

HOST_CORES=$(nproc)

case "$MAX_CORE_PER_JOB" in
    ''|*[!0-9]*) MAX_CORE_PER_JOB=$HOST_CORES ;;
esac
if [ "$MAX_CORE_PER_JOB" -lt 1 ]; then
    MAX_CORE_PER_JOB=$HOST_CORES
fi

CONCURRENCY=$((HOST_CORES / MAX_CORE_PER_JOB))
if [ "$CONCURRENCY" -lt 1 ]; then
    CONCURRENCY=1
fi

export MAX_CORE_PER_JOB

echo "[start_celery_worker] host_cores=$HOST_CORES max_core_per_job=$MAX_CORE_PER_JOB concurrency=$CONCURRENCY"

exec celery -A app_v1.app:celery_app worker -l info --concurrency="$CONCURRENCY"
