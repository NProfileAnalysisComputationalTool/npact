#!/usr/bin/env python
import logging
import os
import os.path
import tempfile
import shutil
import sys
from optparse import OptionParser
from contextlib import contextmanager
import math
from path import path as pathlib
from subprocess import PIPE
from pynpact import parsing
from __init__ import binfile, DATAPATH
import capproc
import prepare
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


def process(verb, filename, config=None, executor=None, outputdir=None):
    if config is None:
        config = parsing.initial(filename, outputdir=outputdir)
    assert executor
    planner = resolve_verb(verb)
    planner(config, executor)
    return config


def run_cmdline(gbkfile):
    sm = None
    try:
        from taskqueue import get_ServerManager
        sm = get_ServerManager(make_server=True, logger=True)
        logging.info("Opening a socket at %s", sm.address)
        sm.start()
        executor = sm.Server()
        config = process('allplots', gbkfile, executor=executor)
        jid = config['pdf_filename']
        logging.info("Work scheduled, waiting")
        executor.get_task(jid).wait()
        logging.info("Finished processing %r", gbkfile)
        logging.info("See output at %r", jid)
    finally:
        if sm:
            sm.shutdown()


if __name__ == '__main__':
    parser = OptionParser("""%prog <genebank file>""")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        exit(1)

    original = args[0]
    gbkfile = os.path.realpath(original)
    logging.basicConfig(
        level=(options.verbose and logging.DEBUG or logging.INFO),
        format="%(asctime)s %(name)-10s %(levelname)-8s %(message)s",
        datefmt='%H:%M:%S')

    run_cmdline(gbkfile)
