import multiprocessing

import random, string


def randomid(length=16):
   return ''.join(random.choice(string.lowercase) for i in range(length))


class InlineExecutor(object):
    "Does everything in line"
    pool = None
    tasks = None

    def __init__(self):
        #self.pool = multiprocessing.Pool()
        self.tasks = {}

    def enqueue(self, callable, id=None, after=None):
        if id is None:
            id = randomid()
        if after is not None:
            for aid in after:
                assert aid in self.tasks, \
                    "The InlineExecutor can't be after a task that doesn't exist yet"

        if id not in self.tasks:
            self.tasks[id] = callable()

        return self.tasks[id]
