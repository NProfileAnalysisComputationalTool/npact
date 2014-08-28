import os.path
import pytest
import py

from taskqueue import get_ServerManager


@pytest.fixture()
def async_executor(request):
    sm = get_ServerManager(make_server=True)
    sm.start()
    request.addfinalizer(sm.shutdown)
    return sm.Server()



class NullExecutor(object):
    "An executor that doens't actually execute anything, just keeps track"
    tasks = None

    def __init__(self):
        self.tasks = {}

    def enqueue(self, callable, tid=None, after=None):
        if tid is None:
            tid = randomid()

        if after is not None:
            for aid in after:
                assert aid in self.tasks, \
                    "The NullExecutor can't be after a task that doesn't exist yet"

        if tid not in self.tasks:
            self.tasks[tid] = callable
        return tid


@pytest.fixture
def null_executor(request):
    return NullExecutor()
