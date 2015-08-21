from __future__ import absolute_import
import mmap
import os
import os.path
import logging
import re
from Bio import SeqIO, Seq
from path import path

from pynpact import genbank
from pynpact.util import mkstemp_rename

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
        'significance': 0.01,

        # allplots
        'first_page_title': None,

        'startBase': 0,
        'basesPerGraph': 10000,
        'graphsPerPage': 5,
        'x-tics': 1000
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
            config['length'] = len(seq)
            ddna = derive_filename(config, filename.getmtime(), 'ddna')
            if not ddna.exists():
                with mkstemp_rename(ddna) as f:
                    f.write(seq.upper())
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
    else:
        config['length'] = int(config['length'])
    return config['length']


def ddna(config):
    if 'ddna' not in config:
        detect_format(config)
    return config['ddna']


def startBase(config):
    if 'startBase' not in config:
        config['startBase'] = 0
    else:
        config['startBase'] = int(config['startBase'])


def endBase(config):
    if 'endBase' not in config:
        config['endBase'] = length(config)
    else:
        config['endBase'] = int(config['endBase'])
    return config['endBase']


def number(val):
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        pass
    return None


def tobool(val):
    if val in (u'false', 0, u'0', False, None):
        return False
    elif val in (u'true', 1, u'1', True):
        return True
    else:
        raise ValueError("Can't convert %r to boolean" % (val,))


def significance(config):
    if 'significance' in config:
        config['significance'] = float(config['significance'])
    else:
        config['significance'] = 0.01
    return config['significance']


def mycoplasma(config):
    if 'mycoplasma' in config:
        config['mycoplasma'] = tobool(config['mycoplasma'])
    else:
        config['mycoplasma'] = False
    return config['mycoplasma']


def first_page_title(config):
    if not config.get('first_page_title'):
        detect_format(config)
        config['first_page_title'] = config.get('description', 'Unknown Genome')
    return config['first_page_title']


def following_page_title(config):
    if 'following_page_title' not in config:
        config['following_page_title'] = 'Page {0}'
    return config['following_page_title']


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
    'startBase': "The base pair coordinate at which to start graphing.",
    'endBase': "The base pair coordinate at which to end graphing.",
    'length': "The length, in base pairs, of the genome being analyzed.",
    'first_page_title': "The title of the page containing the beginning of the genome.",
    'following_page_title': "The title of the pages after the first. Use {0} to get the page number.",
    'skip_prediction': "Should the acgt_gamma prediction be run? (click here to skip)",
    'significance': "What should the acgt_gamma prediction consider significant?",
    'nucleotides': "The bases to count the frequency of on the primary strand.",
    'alternate_colors': "Modifies graph colors from RGB to an alternative set that should be more easily distinguishable by people with certain color-blindness conditions.",
    }


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
#     endBase = endBase


def translate(config, rc=False):
    table = 1
    if mycoplasma(config):
        # table 4 is for mycoplasma ala:
        # http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
        table = 4
    fd, fmap = None, None
    try:
        log.debug("Doing translation with table %d, rc: %s", table, rc)
        fd = os.open(ddna(config), os.O_RDONLY)
        fmap = mmap.mmap(fd, 0, mmap.MAP_SHARED, mmap.PROT_READ)
        # NProfiler.ddna is 0 indexed; the dna by convention
        # (e.g. from the C or NCBI) is 1 indexed.
        startIdx = config['startBase'] - 1
        # The end index here is inclusive but array.slice isn't so we
        # don't need to subtract 1
        endIdx = config['endBase']
        seq = Seq.Seq(fmap[startIdx:endIdx])

        if rc:
            seq = seq.reverse_complement()
        return {
            'seq': str(seq),
            'trans': str(Seq.translate(seq, table))
        }
    finally:
        if fmap:
            fmap.close
        if fd:
            os.close(fd)
