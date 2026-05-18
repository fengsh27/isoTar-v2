"""Per-job and global concurrency limits.

Defaults are tuned for a 32-core CPU budget (4 concurrent jobs x 8 cores/job).
Override at deploy time via env vars:

    ISOTAR_MAX_CORES_PER_JOB     (default 8)   -- API rejects cores > this
    ISOTAR_DEFAULT_CORES         (default 8)   -- used when client omits "cores"
    ISOTAR_MAX_CONCURRENT_JOBS   (default 4)   -- celery --concurrency
"""

import os


def _read_int_env(name, default):
    raw = os.environ.get(name)
    if raw is None or raw == "":
        return default
    try:
        v = int(raw)
        return v if v > 0 else default
    except (TypeError, ValueError):
        return default


MAX_CORES_PER_JOB    = _read_int_env("ISOTAR_MAX_CORES_PER_JOB",    8)
DEFAULT_CORES        = _read_int_env("ISOTAR_DEFAULT_CORES",        8)
MAX_CONCURRENT_JOBS  = _read_int_env("ISOTAR_MAX_CONCURRENT_JOBS",  4)


def validate_cores(value, max_cores=None, default_cores=None):
    """Validate and normalise a client-supplied cores value.

    Returns (cores, error) where error is None on success, or
    a (message, status_code) tuple on failure. When value is None,
    the default is used.

    Accepts ints and integer-valued strings ("8"). Rejects floats,
    non-numeric strings, and out-of-range integers.
    """
    if max_cores is None:
        max_cores = MAX_CORES_PER_JOB
    if default_cores is None:
        default_cores = DEFAULT_CORES

    if value is None:
        return default_cores, None

    # bool is a subclass of int -- reject explicitly so True/False don't pass.
    if isinstance(value, bool):
        return None, ("cores must be a positive integer", 400)

    if isinstance(value, int):
        cores = value
    elif isinstance(value, str):
        try:
            cores = int(value.strip())
        except (TypeError, ValueError):
            return None, ("cores must be a positive integer", 400)
    else:
        return None, ("cores must be a positive integer", 400)

    if cores < 1 or cores > max_cores:
        return None, (
            "cores must be between 1 and {}".format(max_cores),
            400,
        )
    return cores, None
