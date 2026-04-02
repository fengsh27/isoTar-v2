import os

from celery import Celery

BROKER_URL = os.environ.get("CELERY_BROKER_URL", "amqp://")
BACKEND_URL = os.environ.get("CELERY_RESULT_BACKEND", "rpc://")

celery_app = Celery("app_v1", broker=BROKER_URL, backend=BACKEND_URL)
celery_app.conf.update(
    task_serializer="json",
    result_serializer="json",
    accept_content=["json"],
    task_track_started=True,
)
