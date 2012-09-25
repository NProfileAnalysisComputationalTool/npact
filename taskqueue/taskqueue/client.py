import logging
from contextlib import contextmanager
from multiprocessing.managers import RemoteError, SyncManager

import taskqueue


class Manager(SyncManager):
    pass

logger = logging.getLogger(__name__)

Manager.register('the_server')

@contextmanager
def server_call():
    manager = Manager(taskqueue.LISTEN_ADDRESS, authkey=taskqueue.AUTH_KEY)
    manager.connect()
    try:
        yield manager.the_server()
    except RemoteError:
        logger.exception("Error in remote process")
        raise




def ready(id):
    """Check to see if the task with the given id is finished
    executing yet.
    """
    with server_call() as server:
        return server.ready(id)


def enqueue(fn, args=[], kwargs={}):
    """Enqueue the function in the task queue and return an identifier
    that can be used to check status or get the result later.
    """
    with server_call() as server:
        return server.enqueue(fn, args, kwargs)

def result(id):
    """Get the result of the task with the given ID. NB: If the task
    isn't finished this will raise a timeout exception.
    """
    with server_call() as server:
        return server.result(id)

def get_task(id):
    with server_call() as server:
        return server.get_task(id)


def perform(fn, args=[], kwargs={}):
    return get_task(enqueue(fn, args, kwargs)).get()
