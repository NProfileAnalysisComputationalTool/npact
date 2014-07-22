#!/usr/bin/env python
import logging
import os
import os.path
import tempfile
import shutil
import sys
from collections import namedtuple
from optparse import OptionParser
from contextlib import contextmanager
import math
from subprocess import PIPE
from pynpact import parsing
from __init__ import binfile, DATAPATH
import capproc
import prepare
import util

logger = logging.getLogger('pynpact')

# Results = namedtuple('Results', [
#     'File_of_published_accepted_CDSs',
#     'File_list_of_nucleotides_in_200bp_windows',
#     'File_of_new_CDSs',
#     'File_of_published_rejected_CDSs'
#     'File_of_GC_coding_potential_regions',
#     'acgt_gamma_output'
#     'combined_ps_name', 'pdf_filename'])


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
