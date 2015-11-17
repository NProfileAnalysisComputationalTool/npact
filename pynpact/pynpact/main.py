import logging
import os
import os.path
import tempfile
import shutil
import sys
from optparse import OptionParser
from contextlib import contextmanager
from path import Path
from subprocess import PIPE
from pynpact import parsing, executors
from __init__ import binfile, DATAPATH
import capproc
import util

logger = logging.getLogger('pynpact')

from pynpact.steps import extract, allplots, acgt_gamma, nprofile


def resolve_verb(verb):
    "Convert the verb into something that has the plan function"
    mod = {'extract': extract,
           'nprofile': nprofile,
           'allplots': allplots,
           'acgt_gamma': acgt_gamma}[verb]
    return mod.plan


def process(verb, config, executor=None, outputdir=None):
    assert executor
    planner = resolve_verb(verb)
    planner(config, executor)
    return config


@contextmanager
def getexecutor(name):
    sm = None
    name = name.lower()
    try:
        if name == 'server':
            from taskqueue import get_ServerManager
            sm = get_ServerManager(make_server=True, logger=True)
            logging.info("Opening a socket at %s", sm.address)
            sm.start()
            e = sm.Server()
        elif name == 'inline':
            e = executors.InlineExecutor()
        elif name == 'gevent':
            e = executors.GeventExecutor()
        elif name == 'daemon':
            from taskqueue import client
            e = client.get_server()
        else:
            raise ValueError("Unknown executor: %r" % name)
        yield e
    finally:
        if sm:
            sm.shutdown()


def run_cmdline(gbkfile, executorName):
    config = parsing.initial(gbkfile)
    with getexecutor(executorName) as executor:
        config = process('allplots', config, executor=executor)
        jid = config.get('pdf_filename') or config.get('combined_ps_name')
        if not Path(jid).exists():
            logging.info("Work scheduled, waiting")
            output = executor.result(jid, timeout=None)
            logging.info("Finished processing %r", gbkfile)
        else:
            output = jid
        logging.info("See output at %r", output)
