"""Go through the genbank record pulling out gene names and locations

        $ extract MYCGE.gbk 0 gene 0 locus_tag > MYCGE.genes
"""
from __future__ import absolute_import
import logging
import os.path
import sys

from pynpact import binfile, InvalidGBKException
from pynpact import capproc, parsing
from pynpact.util import Hasher, reducedict, mkstemp_rename, delay


logger = logging.getLogger('pynpact.steps.extract')
statuslog = logging.getLogger('pynpact.statuslog')

BIN = binfile("extract")

KEYS = ['GeneDescriptorKey1', 'GeneDescriptorKey2',
        'GeneDescriptorSkip1', 'GeneDescriptorSkip2',
        'filename']

OUTPUTKEY = 'File_of_published_accepted_CDSs'


def plan(config, executor):
    if parsing.isgbk(config):
        rconfig, hash = get_hash(config)
        target_file = parsing.derive_filename(config, hash, 'genes')
        config[OUTPUTKEY] = target_file
        if target_file.exists():
            return None
        executor.enqueue(
            delay(_extract)(config['filename'], target_file, rconfig),
            tid=target_file
        )
        return target_file
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
    if target_file.exists():
        return target_file
    with mkstemp_rename(target_file) as out:
        statuslog.debug(
            "Extracting genes in %s.", os.path.basename(gbkfile))
        cmd = [BIN, gbkfile,
               config['GeneDescriptorSkip1'], config['GeneDescriptorKey1'],
               config['GeneDescriptorSkip2'], config['GeneDescriptorKey2']]
        capproc.capturedCall(
            cmd, check=True, logger=logger, stdout=out, stderr=sys.stderr)
    return target_file
