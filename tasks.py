# function calls to give Celery go here
# for the most part should probably be stuff pulled from urls.py
# after tasks are added, celery can be started from the main thread
# with functions called and sent to the "broker"
# starting with 1 worker (local machine), can add other threads (?) and potentially other machines
from celery import Celery
import os


os.environ.setdefault("DJANGO_SETTINGS_MODULE", "myPhyloDB.settings")

from django.conf import settings

app = Celery('proj')

app.config_from_object('django.conf:settings')
app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)

@app.task
def add(x, y):
    return x + y
