import os
import logging
import tempfile
import sys
from multiprocessing.connection import Client

from optparse import OptionParser

from path import path
import taskqueue


logger = logging.getLogger(__name__)

def client_call(action, **kwargs):
    client = Client(taskqueue.LISTEN_ADDRESS, authkey=taskqueue.AUTH_KEY)
    kwargs['action'] = action
    client.send(kwargs)
    response =  client.recv()
    if response['status'] == 'ok':
        return response['result']
    elif 'exception' in response:
        raise response['exception']
    else:
        raise Exception("Bad Response", response)


def ready(id):
    """Check to see if the task with the given id is finished
    executing yet.
    """
    return client_call('ready', id=id)


def enqueue(fn, args=[], kwargs={}):
    """Enqueue the function in the task queue and return an identifier
    that can be used to check status or get the result later.
    """
    return client_call('enqueue', fn=fn, args=args, kwargs=kwargs)

def result(id):
    """Get the result of the task with the given ID. NB: If the task
    isn't finished this will raise a timeout exception.
    """
    return client_call('result', id=id)
