import os
import logging
import tempfile
import sys
from multiprocessing.connection import Client

from optparse import OptionParser

from path import path
from pynpact import capproc
import taskqueue


logger = logging.getLogger(__name__)


def ready(id):
    client = Client(taskqueue.LISTEN_ADDRESS, authkey=taskqueue.AUTH_KEY)
    client.send({'action': 'ready',
                 'id': id})
    return client.recv()

def enqueue(task):
    client = Client(taskqueue.LISTEN_ADDRESS, authkey=taskqueue.AUTH_KEY)
    client.send({'action': 'enqueue',
                 'task': task})
    return client.recv()
    

def result(id):
    client = Client(taskqueue.LISTEN_ADDRESS, authkey=taskqueue.AUTH_KEY)
    client.send({'action': 'result',
                 'id': id})
    return client.recv()
