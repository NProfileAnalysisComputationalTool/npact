"""Go through the genbank record pulling out gene names and locations

        $ extract MYCGE.gbk 0 gene 0 locus_tag > MYCGE.genes
"""
from __future__ import absolute_import
import logging
import os.path
import sys

from pynpact.steps import derive_filename
from pynpact import binfile
from pynpact import capproc, parsing
from pynpact.util import Hasher, reducedict, mkstemp_overwrite, delay


logger = logging.getLogger('pynpact.steps.extract')
statuslog = logging.getLogger('pynpact.statuslog')

BIN = binfile("extract")

KEYS = ['GeneDescriptorKey1', 'GeneDescriptorKey2',
        'GeneDescriptorSkip1', 'GeneDescriptorSkip2',
        'filename']


def plan(config):
    if parsing.isgbk(config):
        config, hash = get_hash(config)
        target_file = derive_filename(config, hash, 'genes')
        yield (
            delay(_extract)(config['filename'], target_file, config),
            target_file,
            None
        )
    else:
        raise InvalidGBKException()


def get_hash(config):
    config = reducedict(config, KEYS)
    h = Hasher()
    h.hashdict(config)
    h.hashfiletime(BIN)
    h.hashfiletime(config['filename'])
    return config, h.hexdigest()


def _extract(gbkfile, target_file, config):
    with mkstemp_overwrite(target_file) as out:
        statuslog.debug(
            "Extracting genes in %s.", os.path.basename(gbkfile))
        cmd = [BIN, gbkfile,
               config['GeneDescriptorSkip1'], config['GeneDescriptorKey1'],
               config['GeneDescriptorSkip2'], config['GeneDescriptorKey2']]
        capproc.capturedCall(
            cmd, check=True, logger=logger, stdout=out, stderr=sys.stderr)
    return target_file
