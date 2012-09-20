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
    return client_call('ready', id=id)


def enqueue(fn, args=[], kwargs={}):
    return client_call('enqueue', fn=fn, args=args, kwargs=kwargs)

def result(id):
    return client_call('result', id=id)
