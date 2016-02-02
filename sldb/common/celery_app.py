from celery import Celery
from celery.exceptions import TimeoutError

app = Celery('sldb', backend='redis://', broker='amqp://guest@localhost//')

def get_result(task, max_tries=None, timeout=.5):
    tries = 0
    while True:
        if max_tries is not None and tries > max_tries:
            raise TimeoutError()

        try:
            return task.get(timeout=timeout)
        except TimeoutError:
            tries += 1
