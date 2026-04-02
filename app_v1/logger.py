import logging
import os
from logging.handlers import TimedRotatingFileHandler

LOG_DIR   = os.environ.get("ISOTAR_LOG_DIR",   "/opt/out")
LOG_LEVEL = os.environ.get("ISOTAR_LOG_LEVEL", "INFO").upper()

_FORMATTER = logging.Formatter(
    "%(asctime)s %(levelname)-8s [%(name)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def _log_namer(default_name):
    """Rename rotated files from app_v1.log.YYYY-MM-DD → app_v1.YYYY-MM-DD.log."""
    base, date_suffix = default_name.rsplit(".", 1)
    root, ext = base.rsplit(".", 1)
    return "{}.{}.{}".format(root, date_suffix, ext)


def get_logger(name="app_v1"):
    """Return a logger that writes to date-based log files under LOG_DIR.

    Active log:   <LOG_DIR>/app_v1.log
    Rotated logs: <LOG_DIR>/app_v1.YYYY-MM-DD.log  (kept for 30 days)

    Both the Flask process and Celery workers call this; the check on
    logger.handlers prevents duplicate handlers when the module is
    imported multiple times in the same process.
    """
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger

    logger.setLevel(getattr(logging, LOG_LEVEL, logging.INFO))
    logger.propagate = False

    # --- file handler ---
    os.makedirs(LOG_DIR, exist_ok=True)
    log_file = os.path.join(LOG_DIR, "app_v1.log")
    file_handler = TimedRotatingFileHandler(
        log_file,
        when="midnight",
        backupCount=30,
        encoding="utf-8",
    )
    file_handler.suffix = "%Y-%m-%d"
    file_handler.namer = _log_namer
    file_handler.setFormatter(_FORMATTER)
    logger.addHandler(file_handler)

    # --- console handler ---
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(_FORMATTER)
    logger.addHandler(console_handler)

    return logger
