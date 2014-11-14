from __future__ import absolute_import
import os.path
import logging
import re
from Bio import SeqIO
from path import path

from pynpact import genbank

log = logging.getLogger(__name__)

defaults = {
        'nucleotides': ['C', 'G'],

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


def detect_format(config):
    if 'format' not in config:
        try:
            filename = path(config['filename'])
            with filename.open('r') as f:
                header = f.read(5)
            if filename.ext in ('.gb', '.gbk') or header.startswith('LOCUS'):
                log.debug("Attempting %s as genbank", filename)
                seqrec = genbank.parse_seq_rec(config['filename'])
                config['format'] = 'genbank'
                config['id'] = seqrec.id
                config['description'] = seqrec.description
                seq = str(seqrec.seq)
            elif filename.ext in ('.fna', '.fasta') or header.startswith('>'):
                seqrec = SeqIO.read(filename, 'fasta')
                config['format'] = 'fasta'
                config['id'] = seqrec.id
                config['description'] = seqrec.description
                seq = str(seqrec.seq)
            else:
                with filename.open('r') as f:
                    seq = f.read()
                seq = re.sub('\s', '', seq)
                config['format'] = 'raw'
            seq = seq.upper()
            config['length'] = len(seq)
            ddna = derive_filename(config, filename.getmtime(), 'ddna')
            with ddna.open('w') as f:
                f.write(seq)
            config['ddna'] = ddna
        except:
            log.exception("Error detecting format")
            config['format'] = None


def isgbk(config):
    if 'format' not in config:
        detect_format(config)
    return config['format'] == 'genbank'


def isfasta(config):
    if 'format' not in config:
        detect_format(config)
    return config['format'] == 'fasta'


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
        detect_format(config)
        config['first_page_title'] = config.get('description', 'Page 1')
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
