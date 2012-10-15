import logging
from contextlib import contextmanager
import multiprocessing
from multiprocessing.managers import RemoteError, SyncManager



import taskqueue
from taskqueue import tqdaemon

log = logging.getLogger(__name__)


class Manager(SyncManager):
    pass

Manager.register('the_server')


def ensure_daemon():
    """Check to see if daemon is running, attempt to start it if not.

    Timing reported this at ~50 us; shouldn't be bad to keep in main path
    """
    if not tqdaemon.status():
        log.info("Attempting to start missing tqdaemon.")
        multiprocessing.Process(target=tqdaemon.daemonize).start()

@contextmanager
def server_call():
    """Connect to the server and return that object.

    Does a new connection every time. This ensure that if the daemon
    was restarted it still works and that we don't build up any stale
    data. Connection is ~40ms so not too bad.
    """
    ensure_daemon()
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
    "Get the log of the process starting at given position."
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
    "Get a async_result proxy reference to the task"
    with server_call() as server:
        return server.get_task(id)

def perform(fn, args=None, kwargs=None):
    "Do the given task synchronously in the separate daemon process."
    return get_task(enqueue(fn, args, kwargs)).get()

def after(id, fn, additional_args=None, kwargs=None):
    "enqueue new task that will wait on this task and then call the fn"
    return enqueue(_after, [id, fn, additional_args, kwargs])

def _after(id, fn, additional_args=None, kwargs=None):
    task = get_task(id)
    result = task.get()
    if additional_args is None : additional_args = []
    if kwargs is None          : kwargs = {}
    return fn(result, *additional_args, **kwargs)
