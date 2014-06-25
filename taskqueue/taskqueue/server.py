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
import multiprocessing
import os.path
import tempfile
import collections
import sys

from multiprocessing.managers import SyncManager

import taskqueue


log = logging.getLogger(__name__)


class Manager(SyncManager):
    pass


method_to_typeid = {
    'get_task': 'AsyncResult'
}


def start_everything():
    try:

        # it seems that the callable in register is called whenever a
        # new client connects, so need to instantiate it once and
        # always return that in order to have it actually shared.
        the_server = Server()
        Manager.register(
            'the_server',
            callable=lambda: the_server, method_to_typeid=method_to_typeid)
        log.info("Opening a socket at %s", taskqueue.LISTEN_ADDRESS)
        manager = Manager(taskqueue.LISTEN_ADDRESS, authkey=taskqueue.AUTH_KEY)
        manager.get_server().serve_forever()
    except:
        log.exception("Error starting taskqueue server")


class Task(object):
    func = None
    args = None
    kwargs = None

    tid = None
    path = None

    def __repr__(self):
        return "#<Task(func=%r, tid=%r)>" % (self.func, self.tid)

    def __init__(self, func, args=None, kwargs=None, tid=None, path=None):
        self.func = func
        self.args = args
        self.kwargs = kwargs
        self.tid = tid
        self.path = path

    def run(self):
        args = [] if self.args is None else self.args
        kwargs = {} if self.kwargs is None else self.kwargs
        return self.func(*args, **kwargs)

    def mark_finished(self):
        if not self.path:
            log.warning("Can't mark task without path: %r", self.tid)
            return

        try:
            if self.todop:
                os.rename(self.path, self.task_base)
                self.path = self.task_base
            else:
                log.warning("Refinished finished task %r", self.tid)
        except:
            log.exception("Error marking task finished: %r", self.tid)

    @property
    def task_base(self):
        if self.path:
            return os.path.splitext(self.path)[0]

    @property
    def todop(self):
        return self.path and os.path.exists(self.path) and \
            os.path.splitext(self.path)[1] == '.todo'


class Server(object):
    pool = None
    tasks = None
    work_dir = None
    manager = None

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        if not self.work_dir:
            # storage for tasks in progress/log info.
            self.work_dir = os.path.join(taskqueue.BASE_DIR, 'queue')
        if not os.path.exists(self.work_dir):
            os.makedirs(self.work_dir)
        self.tasks = dict()
        self.pool = multiprocessing.Pool()
        log.info("Finished initializing server")

    def ready(self, tid):
        log.debug("Checking on status of %s", tid)
        promise = self.get_task(tid)
        return promise.ready()

    def result(self, tid):
        log.debug("Checking on result of %s", tid)
        promise = self.get_task(tid)
        return promise.get(0.05)

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
        promise = self.pool.apply_async(async_wrapper, [task])
        self.tasks[task.tid] = promise
        return promise

    def enqueue(self, fn, args=None, kwargs=None):
        task = Task(fn, args, kwargs)
        self.pickle_task(task)
        self._enqueue(task)
        return task.tid

    def get_task(self, tid):
        if tid in self.tasks:
            return self.tasks[tid]
        else:
            log.debug("Missing task %s, checking filesystem.", tid)
            task = self.unpickle_task(tid)
            if task.todop:
                return self._enqueue(task)

    def pickle_task(self, task):
        "Writes the task def to a file via pickle in case of server restart"
        with tempfile.NamedTemporaryFile(
                delete=False, dir=self.work_dir,
                prefix='tq-', suffix='.todo') as f:
            log.info("tempfile: %r", f.name)
            cPickle.dump(task, f, cPickle.HIGHEST_PROTOCOL)
            path = f.name
        tid = os.path.splitext(os.path.basename(path))[0]
        task.tid, task.path = tid, path
        return tid, path

    def unpickle_task(self, tid):
        "Read a task definition from file."
        path = os.path.join(self.work_dir, tid)
        if not os.path.exists(path):
            path = path + ".todo"
        if not os.path.exists(path):
            raise taskqueue.NoSuchTaskError()
        with open(path, 'rb') as f:
            task = cPickle.load(f)
        task.tid, task.path = tid, path
        return task


def async_wrapper(task):
    """Wrapper func around executing a task; ensures env setup and teardown.

    * This is run inside a pool process via the pool apply_async above.
    * Ensures the job 'todo' file is marked finished
    * sets up logging to a file for this task
    """

    ### Logging setup:
    ### In this proc all logging should go to a file named after the task

    root_logger = logging.root
    # blow away any other handlers from previous uses of this process.
    root_logger.handlers = []

    log_path = task.task_base + ".log"
    file_handler = logging.FileHandler(log_path, mode='a')
    file_handler.setFormatter(
        logging.Formatter(
            '%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
            datefmt='%Y%m%d %H:%M:%S'))

    root_logger.setLevel(logging.DEBUG)
    root_logger.addHandler(file_handler)

    stderrlog = task.task_base + ".stderr"
    oldstderr = sys.stderr
    with open(stderrlog, 'w') as f:
        sys.stderr = f
        try:
            rc = task.run()
            log.info("Finished %r", task.tid)
        finally:
            task.mark_finished()
        return rc
    sys.stderr = oldstderr
