#!/usr/bin/env python
"""


# ENQUEUE:  mkstemp should guarantee uniquenames

# DEQUEUE: list files in directory, sort by date modified, start working on job.
#  * multiple pool processes?
#   * have one process doing the dequeue, use multiprocessing module beyond that.


http://pypi.python.org/pypi/python-daemon/
"""

import logging
import multiprocessing
import os
import random
import sys
import tempfile

from multiprocessing.connection import Listener
from path import path


import taskqueue



log = logging.getLogger(__name__)

class Task(object):
    backing_file = None
    desc = None
    promise = None

class Server(object):
    socket = None
    pool = None
    tasks = None

    def __init__(self):
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

    def random_id(self):
        id = random.randint(0,16) #TODO:(0,16777216)
        while id in self.tasks:
            id = random.randint(0,16) #TODO:(0,16777216)
        return id

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
                log.exception("Error handling request.")
                return {'status': 'error', 'exception': e}
        else:
            return Exception("Unknown command")

    
    def ready(self, id):
        if id in self.tasks:
            promise = self.tasks[id]
            return promise.ready()
        else:
            raise Exception("No such task")

    def result(self, id):
        if id in self.tasks:
            promise = self.tasks[id]
            return promise.get(0)
        else:
            return Exception("No such task")

    def enqueue(self, task):
        #an AsyncResult object not sure if we can pass it back.
        promise = self.pool.apply_async(*task)
        id = self.random_id()
        self.tasks[id] = promise
        return id
