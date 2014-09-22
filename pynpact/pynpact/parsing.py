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
    return config['first_page_title']



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



CONFIG_HELP_TEXT = {
    'start_base': "The base pair coordinate at which to start graphing.",
    'end_base': "The base pair coordinate at which to end graphing.",
    'length': "The length, in base pairs, of the genome being analyzed.",
    'first_page_title': "The title of the page containing the beginning of the genome.",
    'following_page_title': "The title of the pages after the first. Use {0} to get the page number.",
    'skip_prediction': "Should the acgt_gamma prediction be run? (click here to skip)",
    'significance': "What should the acgt_gamma prediction consider significant?",
    'nucleotides': "The bases to count the frequency of on the primary strand.",
    'alternate_colors': "Modifies graph colors from RGB to an alternative set that should be more easily distinguishable by people with certain color-blindness conditions.",
    }

significance_levels = ("0.01", "0.001", "0.0001")
significance_levels = zip(significance_levels, significance_levels)




## playing around with this idea, haven't really gotten it yet
# class Config(object):
#     def __init__(self, filename):
#         self.__dict__ = initial(filename)
#     def __getitem__(self, key):
#         if key in self.__dict__:
#             return self.__dict__[key]
#         elif key in defaults:
#             return defaults[key]
#         elif hasattr(self, key):
#             val = getattr(self, key)()
#             self[key] = val
#             return val
#         else:
#             raise KeyError(key)
#     def __setitem__(self, key, val):
#         self.__dict__[key] = val
#     def __contains__(self, key):
#         return key in self.__dict__
#     first_page_title = first_page_title
#     end_base = end_base
