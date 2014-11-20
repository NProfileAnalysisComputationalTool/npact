"""Go through the genbank record pulling out gene names and locations

        $ extract MYCGE.gbk 0 gene 0 locus_tag > MYCGE.genes
"""
from __future__ import absolute_import
import logging
import os.path
import sys

from pynpact import binfile, InvalidGBKException
from pynpact import capproc, parsing
from pynpact.util import Hasher, reducedict
from pynpact.steps import producer, enqueue


logger = logging.getLogger('pynpact.steps.extract')
statuslog = logging.getLogger('pynpact.statuslog')

BIN = binfile("extract")

KEYS = ['GeneDescriptorKey1', 'GeneDescriptorKey2',
        'GeneDescriptorSkip1', 'GeneDescriptorSkip2',
        'filename']

OUTPUTKEY = 'File_of_published_accepted_CDSs'


def plan(config, executor):
    if parsing.isgbk(config):
        logger.debug(
            "GBK file, extracting known gene names %s", config['filename'])
        rconfig, hash = get_hash(config)
        target_file = parsing.derive_filename(config, hash, 'genes')
        config[OUTPUTKEY] = target_file
        return enqueue(_extract, executor, rconfig, target_file)


def get_hash(config):
    config = reducedict(config, KEYS)
    h = Hasher()
    h.hashdict(config)
    h.hashfiletime(BIN)
    h.hashfiletime(config['filename'])
    return config, h.hexdigest()


@producer()
def _extract(config, out):
    gbkfile = config['filename']
    statuslog.debug(
        "Extracting genes in %s.", os.path.basename(gbkfile))
    cmd = [BIN, gbkfile,
           config['GeneDescriptorSkip1'], config['GeneDescriptorKey1'],
           config['GeneDescriptorSkip2'], config['GeneDescriptorKey2']]
    capproc.capturedCall(
        cmd, check=True, logger=logger, stdout=out, stderr=sys.stderr)
