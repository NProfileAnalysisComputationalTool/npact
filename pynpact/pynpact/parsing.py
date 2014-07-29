from __future__ import absolute_import
import os.path
import logging
from Bio import SeqIO
from path import path

from pynpact import genbank

log = logging.getLogger(__name__)

defaults = {
        'nucleotides': ['c', 'g'],

        # keys for extract.c
        'GeneDescriptorKey1': 'gene',
        'GeneDescriptorKey2': 'locus_tag',
        'GeneDescriptorSkip1': 0,
        'GeneDescriptorSkip2': 0,

        # keys for CG:
        'window_size': 201,
        'step': 51,
        'period': 3,

        # acgt_gamma:
        'skip_prediction': False,
        'significance': "0.001",

        ##allplots
        #http://docs.python.org/library/string.html#formatspec
        'first_page_title': None,
        'following_page_title': 'Page {0}',

        'bp_per_page': 50000,
        'start_base': 0,

        }


def initial(filename, outputdir=None):
    config = defaults.copy()
    assert os.path.exists(filename)
    config['filename'] = filename
    config['basename'] = os.path.basename(filename)
    if outputdir is None:
        outputdir = os.path.dirname(os.path.realpath(filename))
    config['outputdir'] = outputdir
    return config


def isgbk(config):
    if 'isgbk' not in config:
        try:
            seqrec = genbank.parse_seq_rec(config['filename'])
            if not config.get('first_page_title'):
                config['first_page_title'] = seqrec.description
            config['id'] = seqrec.id
            config['length'] = len(seqrec)
            config['isgbk'] = True
        except Exception as e:
            log.debug('Error parsing gbk file: %s', e)
            config['isgbk'] = False
    return config['isgbk']


def isfasta(config):
    if 'isfasta' not in config:
        try:
            seqrec = SeqIO.read(config['filename'], 'fasta')
            config['length'] = len(seqrec)
            config['isfasta'] = True
        except:
            config['isfasta'] = False
    return config['isfasta']


def length(config):
    if 'length' not in config:
        if isgbk(config):
            # isgbk fills length if it is a gbk file
            pass
        elif isfasta(config):
            # isfasta fills length
            pass
        else:
            config['length'] = 0
    return config['length']


def end_base(config):
    if 'end_base' not in config:
        config['end_base'] = length(config)
    return config['end_base']


def first_page_title(config):
    if not config.get('first_page_title'):
        if isgbk(config):
            pass
        else:
            'Page 1'


def derive_filename(config, hash, newext):
    "Build target filename based on identifying pieces"
    if newext[0] == '.':
        newext = newext[1:]
    filename = path(config['filename'])
    if 'outputdir' in config:
        outputdir = path(config['outputdir'])
    else:
        outputdir = filename.dirname()
    newfilename = '%s-%s.%s' % (filename.namebase, hash, newext)
    return outputdir / newfilename
