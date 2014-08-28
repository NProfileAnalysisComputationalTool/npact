import logging
from contextlib import contextmanager
import multiprocessing
from multiprocessing.managers import RemoteError
import time

from taskqueue import get_ServerManager
from taskqueue import tqdaemon

log = logging.getLogger(__name__)

ENSURE_DAEMON = False


def ensure_daemon():
    """Check to see if daemon is running, attempt to start it if not.

    Timing reported this at ~50 us if the daemon was already running;
    shouldn't be bad to keep in main path. If the daemon wasn't
    running waiting for it to start is a good thing.

    """
    if not tqdaemon.status():
        log.info("Attempting to start missing tqdaemon.")
        p = multiprocessing.Process(target=tqdaemon.daemonize)
        p.start()
        p.join(1)
        time.sleep(0.1)


def get_server():
    if ENSURE_DAEMON:
        ensure_daemon()

    manager = get_ServerManager()
    manager.connect()
    return manager.Server()


@contextmanager
def server_call():
    """Connect to the server and return that object.

    Does a new connection every time. This ensure that if the daemon
    was restarted it still works and that we don't build up any stale
    data. Connection is ~40ms so not too bad.
    """
    try:
        s = get_server()
    except:
        log.exception("Error connecting to server")
        raise
    try:
        yield s
    except RemoteError:
        log.exception("Error in remote process")
        raise


def ready(tid):
    """Check to see if the task with the given id is finished
    executing yet.
    """
    with server_call() as server:
        return server.ready(tid)


def log_output(tid, position=0):
    "Get the log of the process starting at given position."
    with server_call() as server:
        return server.log_output(tid, position)


def log_tail(tid, line_count=1):
    with server_call() as server:
        return server.log_tail(tid, line_count)


def enqueue(fn, tid=None, after=None):
    """Enqueue the function in the task queue and return an identifier
    that can be used to check status or get the result later.
    """
    with server_call() as server:
        return server.enqueue(fn, tid=tid, after=after)


def result(tid):
    """Get the result of the task with the given ID. NB: If the task
    isn't finished this will raise a timeout exception.
    """
    with server_call() as server:
        return server.result(tid)


def get_task(tid):
    "Get a async_result proxy reference to the task"
    with server_call() as server:
        return server.get_task(tid)


def perform(fn):
    "Do the given task synchronously in the separate daemon process."
    return get_task(enqueue(fn)).get()


# def after(tid, fn, additional_args=None, kwargs=None):
#     "enqueue new task that will wait on this task and then call the fn"
#     return enqueue(_after, [tid, fn, additional_args, kwargs])


# def _after(tid, fn, additional_args=None, kwargs=None):
#     "helper function for after"
#     task = get_task(tid)
#     task_result = task.get()
#     if additional_args is None:
#         additional_args = []
#     if kwargs is None:
#         kwargs = {}
#     return fn(task_result, *additional_args, **kwargs)
