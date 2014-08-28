import string
import random


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
