#!/usr/bin/env python
"""
This module listens on a TCP socket and runs a worker pool.

The socket adress is specified by ``taskqueue.LISTEN_ADDRESS``.

Tasks (a python callable and arguments) can be enqueued into the pool
which returns an alphanumeric id to calling connection. That id can
then be used to request the current status or final result of the
task.

"""
import logging
import multiprocessing
import os.path
import collections
import sys

from path import path

import taskqueue
from taskqueue import NoSuchTaskError
from taskqueue.task import Task, async_wrapper


log = logging.getLogger(__name__)


class Server(object):
    pool = None
    work_dir = None
    manager = None

    # TODO: this dict needs to become a weakref.WeakValueDict + deque
    # so we don't hang on to every task ever.
    tasks = None

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        if not self.work_dir:
            # storage for tasks in progress/log info.
            self.work_dir = os.path.join(taskqueue.BASE_DIR, 'queue')
        if not os.path.exists(self.work_dir):
            log.debug("Making work_dir: %r", self.work_dir)
            os.makedirs(self.work_dir)
        self.tasks = dict()
        if self.pool is None:
            # can use pool=False to disable
            self.pool = multiprocessing.Pool()
        log.info("Finished initializing server, work_dir: %s", self.work_dir)

    def get_task(self, tid):
        if tid not in self.tasks:
            log.debug("Missing task %s, checking filesystem.", tid)
            task = Task.unpickle(tid, self.work_dir)
            self._enqueue(task)
        return self.tasks[tid]

    def wait(self, tid, timeout=None):
        return self.get_task(tid).wait(timeout)

    def task_count(self):
        return len(self.tasks)

    def ready(self, tid):
        log.debug("Checking on status of %s", tid)
        promise = self.get_task(tid)
        return promise.ready()

    def result(self, tid, timeout=0.01):
        log.debug("Checking on result of %s", tid)
        promise = self.get_task(tid)
        return promise.get(timeout)

    ## TODO: these need to use sanitize_id or a stored log path for the task
    # def log_output(self, tid, position=0):
    #     log.debug("Retrieving log output for %s from pos:%s", tid, position)
    #     pth = os.path.join(self.work_dir, tid + '.log')
    #     with open(pth, 'rb') as f:
    #         f.seek(position)
    #         buf = f.read()
    #     log.debug("Read %d bytes from log.", len(buf))
    #     return buf

    # def log_tail(self, tid, line_count=1):
    #     log.debug("Retrieving %d log lines for %s", line_count, tid)

    #     # make sure we have the running task (raises exception if missing)
    #     self.get_task(tid)

    #     pth = os.path.join(self.work_dir, tid + '.stderr')
    #     if not os.path.exists(pth):
    #         return ["Queued\n"]
    #     try:
    #         with open(pth, 'rb') as f:
    #             try:
    #                 # 2 is from end
    #                 f.seek(- (line_count + 1) * 120, 2)
    #             except:
    #                 f.seek(0)
    #             lines = collections.deque(maxlen=line_count)
    #             while True:
    #                 line = f.readline()
    #                 if line == '':
    #                     return list(lines)
    #                 elif line != '\n':
    #                     lines.append(line)
    #     except:
    #         log.exception("Error reading lines from %r", pth)
    #         raise

    def _enqueue(self, task):
        log.info('Enqueuing %r', task)
        assert task.tid
        promise = self.pool.apply_async(async_wrapper, [task, self.work_dir])
        self.tasks[task.tid] = promise
        return promise

    def enqueue(self, fn, tid=None, after=None):
        if tid:
            if tid in self.tasks:
                return tid
        if after:
            after = filter(None, after)
            # check that all tasks specified exist or raise the
            # exception now if one of them doesn't. an effect this
            # enforces is that everything must be queued into the pool
            # before this task is to prevent deadlock
            map(self.get_task, after)
        task = Task(fn, tid=tid, after=after)
        task.pickle(self.work_dir)
        self._enqueue(task)
        return task.tid
