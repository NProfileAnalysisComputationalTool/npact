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

from pynpact import capproc
import taskqueue



log = logging.getLogger(__name__)


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
                response = {'status': 'error', 'description': e}

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
            raise
        else:
            action = message['action']
        if 'enqueue' == action:
            task = message['task']
            #an AsyncResult object not sure if we can pass it back.
            promise = self.pool.apply_async(*task)
            id = self.random_id()
            self.tasks[id] = promise
            return {'status': 'ok', 'id': id, }

        elif 'ready' == action:
            id = message['id']
            if id in self.tasks:
                promise = self.tasks[id]
                return {'status': 'ok', 'ready': promise.ready(), 'id': id}
            else:
                return Exception("No such task")
        elif 'result' == action:
            id = message['id']
            if id in self.tasks:
                promise = self.tasks[id]
                try:
                    return {'status':'ok',
                            'ready': True,
                            'result': promise.get()}
                except multiprocessing.TimeoutError, te:
                    return {'status': 'ok', 'ready': False, 'id': id}
            else:
                return Exception("No such task")
        else:
            return Exception("Unknown command")
