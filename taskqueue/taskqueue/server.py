#!/usr/bin/env python
"""
This module listens on a TCP socket and runs a worker pool.

The socket adress is specified by ``taskqueue.LISTEN_ADDRESS``.

Tasks (a python callable and arguments) can be enqueued into the pool
which returns an alphanumeric id to calling connection. That id can
then be used to request the current status or final result of the
task.

"""
import cPickle
import logging
import errno
import multiprocessing
import os.path
import collections
import sys
import string
import random
from multiprocessing.managers import SyncManager

from path import path

import taskqueue


log = logging.getLogger(__name__)


def sanitize_id(tid):
    return tid.replace(os.path.sep, '==')


def randomid():
    length = 16
    return ''.join(
        random.choice(string.lowercase) for i in range(length))


class Task(object):
    func = None
    tid = None
    path = None

    def __repr__(self):
        return "#<Task(func=%r, tid=%r)>" % (self.func, self.tid)

    def __init__(self, func, tid=None):
        self.func = func
        self.tid = tid

    def run(self):
        return self.func()

    # def mark_finished(self):
    #     if not self.path:
    #         log.warning("Can't mark task without path: %r", self.tid)
    #         return

    #     try:
    #         if self.unfinished:
    #             os.rename(self.path, self.task_base)
    #             self.path = self.task_base
    #         else:
    #             log.warning("Refinished finished task %r", self.tid)
    #     except:
    #         log.exception("Error marking task finished: %r", self.tid)

    def pickle(self, work_dir):
        "Writes the task def to a file via pickle in case of server restart"

        work_dir = path(work_dir)
        while True:
            fd = None
            # Unless given a specific tid generate one randomly
            tid = self.tid or randomid()
            p = work_dir / sanitize_id(tid)
            log.debug("Pickling %r to %r", self, p)
            try:
                # write out the file, make sure we have unique handle to it
                fd = os.open(p, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                fd = os.fdopen(fd, 'w')
                cPickle.dump(self.func, fd, cPickle.HIGHEST_PROTOCOL)
                self.tid = tid
            except OSError as e:
                if e.errno == errno.EEXIST:
                    # file wasn't unique; if it was a specified id
                    # then we're already done, if was a random id try
                    # again
                    if tid is self.tid:
                        return self.tid
                    else:
                        log.debug("id collision; retrying")
                        continue
                else:
                    raise
            finally:
                if fd:
                    fd.close()

        return self.tid

    @classmethod
    def unpickle(klass, tid, work_dir):
        "Read a task definition from file."
        path = os.path.join(work_dir, sanitize_id(tid))
        if not os.path.exists(path):
            raise taskqueue.NoSuchTaskError()
        with open(path, 'r') as f:
            func = cPickle.load(f)
        return klass(func=func, tid=tid)


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
        self.pool = multiprocessing.Pool()
        log.info("Finished initializing server")

    def get_task(self, tid):
        if tid in self.tasks:
            return self.tasks[tid]
        else:
            log.debug("Missing task %s, checking filesystem.", tid)
            task = Task.unpickle(tid, self.work_dir)
            if task.unfinished:
                return self._enqueue(task)

    def ready(self, tid):
        log.debug("Checking on status of %s", tid)
        promise = self.get_task(tid)
        return promise.ready()

    def result(self, tid):
        log.debug("Checking on result of %s", tid)
        promise = self.get_task(tid)
        return promise.get(0.01)

    def log_output(self, tid, position=0):
        log.debug("Retrieving log output for %s from pos:%s", tid, position)
        path = os.path.join(self.work_dir, tid + '.log')
        with open(path, 'rb') as f:
            f.seek(position)
            buf = f.read()
        log.debug("Read %d bytes from log.", len(buf))
        return buf

    def log_tail(self, tid, line_count=1):
        log.debug("Retrieving %d log lines for %s", line_count, tid)

        # make sure we have the running task (raises exception if missing)
        self.get_task(tid)

        path = os.path.join(self.work_dir, tid + '.stderr')
        if not os.path.exists(path):
            return ["Queued\n"]
        try:
            with open(path, 'rb') as f:
                try:
                    # 2 is from end
                    f.seek(- (line_count + 1) * 120, 2)
                except:
                    f.seek(0)
                lines = collections.deque(maxlen=line_count)
                while True:
                    line = f.readline()
                    if line == '':
                        return list(lines)
                    elif line != '\n':
                        lines.append(line)
        except:
            log.exception("Error reading lines from %r", path)
            raise

    def _enqueue(self, task):
        log.info('Enqueuing %r', task)
        assert task.tid
        if task.tid not in self.tasks:
            promise = self.pool.apply_async(
                async_wrapper, [task, self.work_dir])
            self.tasks[task.tid] = promise
        else:
            promise = self.tasks[task.tid]
        return promise

    def enqueue(self, fn, tid=None, after=None):
        task = Task(fn, tid=tid)
        if after:
            map(self.get_task, after)
        task.pickle(self.work_dir)

        self._enqueue(task)
        return task.tid


def async_wrapper(task, work_dir):
    """Wrapper func around executing a task; ensures env setup and teardown.

    * This is run inside a pool process via the pool apply_async above.
    * Ensures the job 'todo' file is marked finished
    * sets up logging to a file for this task
    """
    work_dir = path(work_dir)

    # Logging setup:
    # In this proc all logging should go to a file named after the task
    root_logger = logging.root
    # blow away any other handlers from previous uses of this process.
    root_logger.handlers = []
    log_path = work_dir / (task.tid + '.log')
    file_handler = logging.FileHandler(log_path, mode='a')
    file_handler.setFormatter(
        logging.Formatter(
            '%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
            datefmt='%Y%m%d %H:%M:%S'))

    root_logger.setLevel(logging.DEBUG)
    root_logger.addHandler(file_handler)

    stderrlog = work_dir / (task.tid + '.stderr')
    oldstderr = sys.stderr
    try:
        with open(stderrlog, 'w') as f:
            sys.stderr = f
            rc = task.run()
            log.info("Finished %r", task.tid)
            return rc
    finally:
        sys.stderr = oldstderr


class ServerManager(SyncManager):
    pass

method_to_typeid = {
    'get_task': 'AsyncResult'
}


def start_ServerManager():
    try:
        # it seems that the callable in register is called whenever a
        # new client connects, so need to instantiate it once and
        # always return that in order to have it actually shared.
        the_server = Server()
        ServerManager.register(
            'the_server',
            callable=lambda: the_server, method_to_typeid=method_to_typeid)
        log.info("Opening a socket at %s", taskqueue.LISTEN_ADDRESS)
        manager = ServerManager(
            taskqueue.LISTEN_ADDRESS, authkey=taskqueue.AUTH_KEY)
        manager.get_server().serve_forever()
    except:
        log.exception("Error starting taskqueue server")
