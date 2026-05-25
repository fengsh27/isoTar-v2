"""Resolve the running application version.

Resolution order:
  1. ISOTAR_VERSION env var (set from the build arg in the Docker image).
  2. A bundled VERSION file (copied to /app_v1/VERSION in the image).
  3. The repo-root VERSION file (when running from a source checkout).
  4. "unknown".
"""
import os

_ENV = "ISOTAR_VERSION"
_HERE = os.path.dirname(os.path.abspath(__file__))
_CANDIDATES = (
    os.path.join(_HERE, "VERSION"),            # bundled next to this module (image)
    os.path.join(_HERE, "..", "VERSION"),      # repo root (source checkout)
)


def get_version():
    env = os.environ.get(_ENV)
    if env and env.strip():
        return env.strip()
    for path in _CANDIDATES:
        try:
            with open(path) as f:
                value = f.read().strip()
            if value:
                return value
        except (IOError, OSError):
            continue
    return "unknown"
