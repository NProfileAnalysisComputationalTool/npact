import logging
from contextlib import contextmanager
from multiprocessing.managers import RemoteError, SyncManager

import taskqueue


class Manager(SyncManager):
    pass

log = logging.getLogger(__name__)

Manager.register('the_server')

@contextmanager
def server_call():
    manager = Manager(taskqueue.LISTEN_ADDRESS, authkey=taskqueue.AUTH_KEY)
    manager.connect()
    try:
        yield manager.the_server()
    except RemoteError:
        log.exception("Error in remote process")
        raise




def ready(id):
    """Check to see if the task with the given id is finished
    executing yet.
    """
    with server_call() as server:
        return server.ready(id)

def log_output(id, position=0):
    with server_call() as server:
        return server.log_output(id, position)

def enqueue(fn, args=None, kwargs=None):
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


def perform(fn, args=None, kwargs=None):
    return get_task(enqueue(fn, args, kwargs)).get()

def after(id, fn, additional_args=None, kwargs=None):

    #enqueue new task that will wait on this task and then call the fn
    return enqueue(_after, [id, fn, additional_args, kwargs])
    
def _after(id, fn, additional_args=None, kwargs=None):
    task = get_task(id)
    result = task.get()
    if additional_args is None : additional_args = []
    if kwargs is None          : kwargs = {}
    return fn(result, *additional_args, **kwargs)
