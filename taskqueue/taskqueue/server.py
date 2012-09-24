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
from multiprocessing.connection import Listener

import taskqueue


log = logging.getLogger(__name__)


class Server(object):
    socket = None
    pool = None
    tasks = None
    helper = None
    helper_queue = None
    work_dir = '/tmp/'

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        self.pool = multiprocessing.Pool()
        self.tasks = dict()
        log.info("Opening a socket at %r", taskqueue.LISTEN_ADDRESS)
        self.socket = Listener(taskqueue.LISTEN_ADDRESS, authkey=taskqueue.AUTH_KEY, backlog=4)


    def run(self):
        log.info("Accepting connections")
        while True:
            pipe = self.socket.accept()
            message = pipe.recv()
            try:
                response = self.dispatch(message)
            except Exception,e:
                response = {'status': 'error', 'exception': e}
            pipe.send(response)
        log.info("Exiting")

    def dispatch(self, message):
        log.debug("Got message: %r", message)
        if not 'action' in message:
            raise ArgumentError('Message must contain action')
        else:
            action = message.pop('action')

        if action in ['enqueue','result','ready']:
            try:
                r = getattr(self, action)(**message)
                return {'status': 'ok', 'result': r}
            except Exception,e:
                if not isinstance(e, multiprocessing.TimeoutError):
                    log.exception("Error handling request.")
                return {'status': 'error', 'exception': e}
        else:
            return Exception("Unknown command")


    def ready(self, id):
        promise = self.get_task(id)
        return promise.ready()


    def result(self, id):
        promise = self.get_task(id)
        return promise.get(0.05)

    def enqueue(self, fn, args=[], kwargs={}):
        path = self.pickle_task([fn, args, kwargs])
        id = os.path.splitext(os.path.basename(path))[0]
        self.tasks[id] = self.pool.apply_async(async_wrapper, [id, path, fn, args, kwargs])

        return id

    def get_task(self, id):
        if id in self.tasks:
            return self.tasks[id]
        else:
            log.debug("Missing task %r, checking filesystem.", id)
            path = os.path.join(self.work_dir, id)
            if not os.path.exists(path):
                path = path + ".todo"
            if os.path.exists(path):
                return self.unpickle_task(id, path)
            else:
                raise taskqueue.NoSuchTaskError()

    def pickle_task(self, task):
        "Writes the task definition to a file via pickle in case of server restart"
        with tempfile.NamedTemporaryFile(delete=False, dir=self.work_dir, prefix='tq-', suffix='.todo') as f:
            log.info("tempfile: %r", f.name)
            cPickle.dump(task, f, cPickle.HIGHEST_PROTOCOL)
            return f.name

    def unpickle_task(self, id, path):
        "Read a task definition from file and restart that process."
        with open(path, 'rb') as f:
            task = cPickle.load(f)
        if path.endswith('todo'):
            log.warning("Unpickling a TODO task: %r", path)
        promise = self.pool.apply_async(async_wrapper, [id, path] + task)
        self.tasks[id] = promise
        return promise



def async_wrapper(id, task_path, fn, args, kwargs):
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
    #remove any other handlers from previous uses of this process.
    for h in root_logger.handlers: root_logger.removeHandler(h)
    root_logger.setLevel(logging.DEBUG)
    root_logger.addHandler(file_handler)

    try:
        rc = fn(*args, **kwargs)
        log.debug("Finished %r", id)
    finally:
        try:
            if os.path.exists(task_path) and task_ext == '.todo':
                os.rename(task_path, task_base)
            else:
                log.warning("Refinished finished task %r", id)
        except:
            log.exception("Error cleaning up after finished task: %r", id)
    return rc
