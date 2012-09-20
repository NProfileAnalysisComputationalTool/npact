#!/usr/bin/env python
"""


# ENQUEUE:  mkstemp should guarantee uniquenames

# DEQUEUE: list files in directory, sort by date modified, start working on job.
#  * multiple pool processes?
#   * have one process doing the dequeue, use multiprocessing module beyond that.


http://pypi.python.org/pypi/python-daemon/
"""
import cPickle
import logging
import multiprocessing
import os
import os.path
import random
import sys
import tempfile

from multiprocessing.connection import Listener
from path import path


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


    def create_pool(self):
        pass

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
        self.tasks[id] = self.pool.apply_async(async_do_and_cleanup, [id, path, fn, args, kwargs])

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
                raise Exception("No such task")

    def pickle_task(self, task):
        with tempfile.NamedTemporaryFile(delete=False, dir=self.work_dir, prefix='tq-', suffix='.todo') as f:
            log.info("tempfile: %r", f.name)
            cPickle.dump(task, f, cPickle.HIGHEST_PROTOCOL)
            return f.name

    def unpickle_task(self, id, path):
        with open(path, 'rb') as f:
            task = cPickle.load(f)
        if path.endswith('todo'):
            log.warning("Unpickling a TODO task: %r", path)
        promise = self.pool.apply_async(async_do_and_cleanup, [id, path] + task)
        self.tasks[id] = promise
        return promise

def async_do_and_cleanup(id, task_path, fn, args, kwargs):
    try:
        rc = fn(*args, **kwargs)
        log.debug("Finished %r", id)
    finally:
        try:
            root,ext = os.path.splitext(task_path)
            if os.path.exists(task_path) and ext == '.todo':
                os.rename(task_path, root)
            else:
                log.warning("Refinished finished task %r", id)
        except:
            log.exception("Error cleaning up after finished task: %r", id)
    return rc
