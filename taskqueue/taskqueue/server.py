#!/usr/bin/env python
"""
This module listens on a TCP socket (specified by taskqueue.LISTEN_ADDRESS) and runs a worker pool.

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
        the_server = Server()
        Manager.register('the_server', callable=lambda: the_server, method_to_typeid=method_to_typeid)
        log.info("Opening a socket at %s", taskqueue.LISTEN_ADDRESS)
        manager = Manager(taskqueue.LISTEN_ADDRESS, authkey=taskqueue.AUTH_KEY)
        manager.get_server().serve_forever()
    except:
        log.exception("Error starting taskqueue server")


class Server(object):
    pool = None
    tasks = None
    work_dir = None
    manager = None
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        if not self.work_dir:
            #storage for tasks in progress/log info.
            self.work_dir = os.path.join(taskqueue.BASE_DIR, 'queue')
        if not os.path.exists(self.work_dir):
            os.makedirs(self.work_dir)
        self.tasks = dict()
        self.pool = multiprocessing.Pool()


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
        path = os.path.join(self.work_dir, tid + '.stderr')
        with open(path, 'rb') as f:
            try:
                f.seek(- (line_count + 1) * 120, 2) # 2 is from end
            except:
                f.seek(0)
            lines = collections.deque(maxlen=line_count)
            while True:
                line = f.readline()
                if line == '':
                    return list(lines)
                elif line != '\n':
                    lines.append(line)

    def _enqueue(self, tid, path, task):
        log.info('Enqueuing [%r,...]', task[0])
        promise = self.pool.apply_async(async_wrapper, [tid, path] + task)
        self.tasks[tid] = promise
        return promise

    def enqueue(self, fn, args=None, kwargs=None):
        path = self.pickle_task([fn, args, kwargs])
        tid = os.path.splitext(os.path.basename(path))[0]
        self._enqueue(tid, path, [fn, args, kwargs])

        return tid

    def get_task(self, tid):
        if tid in self.tasks:
            return self.tasks[tid]
        else:
            log.debug("Missing task %s, checking filesystem.", tid)
            path = os.path.join(self.work_dir, tid)
            if not os.path.exists(path):
                path = path + ".todo"
            if os.path.exists(path):
                return self.unpickle_task(tid, path)
            else:
                raise taskqueue.NoSuchTaskError()

    def pickle_task(self, task):
        "Writes the task definition to a file via pickle in case of server restart"
        with tempfile.NamedTemporaryFile(delete=False, dir=self.work_dir, prefix='tq-', suffix='.todo') as f:
            log.info("tempfile: %r", f.name)
            cPickle.dump(task, f, cPickle.HIGHEST_PROTOCOL)
            return f.name

    def unpickle_task(self, tid, path):
        "Read a task definition from file and restart that process."
        with open(path, 'rb') as f:
            task = cPickle.load(f)
        if path.endswith('todo'):
            log.warning("Unpickling a TODO task: %r", path)
        return self._enqueue(tid, path, task)



def async_wrapper(tid, task_path, fn, args, kwargs):
    """Wrapper function around executing a task to ensure env setup and teardown.

    * This is run inside a pool process via the pool apply_async above.
    * Ensures the job 'todo' file is marked finished
    * sets up logging to a file for this task
    """
    task_base,task_ext = os.path.splitext(task_path)

    ### Logging setup:
    ### In this pro all logging should go to a file named after the task

    log_path = task_base + ".log"
    file_handler = logging.FileHandler(log_path, mode='a')
    file_handler.setFormatter(logging.Formatter('%(asctime)s %(name)-10s'
                                                ' %(levelname)-8s %(message)s',
                                                datefmt='%Y%m%d %H:%M:%S'))
    root_logger = logging.getLogger('')
    #blow away any other handlers from previous uses of this process.
    root_logger.handlers = []
    root_logger.setLevel(logging.DEBUG)
    root_logger.addHandler(file_handler)

    stderrlog = task_base + ".stderr"
    oldstderr = sys.stderr
    with open(stderrlog, 'w') as f:
        sys.stderr = f
        try:
            if args is None:   args = []
            if kwargs is None: kwargs = {}
            rc = fn(*args, **kwargs)
            log.debug("Finished %r", tid)
        finally:
            try:
                if os.path.exists(task_path) and task_ext == '.todo':
                    os.rename(task_path, task_base)
                else:
                    log.warning("Refinished finished task %r", tid)
            except:
                log.exception("Error cleaning up after finished task: %r", tid)
        return rc
    sys.stderr = oldstderr
