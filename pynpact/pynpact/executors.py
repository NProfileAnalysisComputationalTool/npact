import string
import random

import logging

log = logging.getLogger(__name__)


def randomid(length=16):
    return ''.join(
        random.choice(string.lowercase) for i in range(length))


class InlineExecutor(object):
    "Does everything in line"
    pool = None
    tasks = None

    def __init__(self):
        self.tasks = {}

    def enqueue(self, callable, tid=None, after=None):
        if tid is None:
            tid = randomid()

        if after is not None:
            for aid in after:
                if aid:
                    assert aid in self.tasks, \
                        "After failed: %s, %s" % (tid, aid)

        if tid not in self.tasks:
            self.tasks[tid] = callable()
        return tid

    def result(self, tid, **kwargs):
        return self.tasks[tid]


class GeventExecutor(object):
    "Does everything in line"
    pool = None
    tasks = None

    def __init__(self):
        from gevent.monkey import patch_all
        patch_all(subprocess=True)
        self.tasks = {}

    def get_task(self, tid):
        return self.tasks[tid]

    def enqueue(self, callable, tid=None, after=None):
        import gevent
        if tid is None:
            tid = randomid()
        if tid in self.tasks:
            return tid

        if after is not None:
            gevent.wait(map(self.get_task, after))

        self.tasks[tid] = gevent.spawn(callable)
        return tid

    def ready(self, tid):
        self.get_task(tid).ready()

    def result(self, tid, timeout=0.01):
        return self.get_task(tid).get(timeout=timeout)
