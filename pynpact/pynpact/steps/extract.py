"""Go through the genbank record pulling out gene names and locations

        $ extract MYCGE.gbk 0 gene 0 locus_tag > MYCGE.genes
"""
from __future__ import absolute_import
import logging
import os.path
import sys

from pynpact.steps import BaseStep, delay
from pynpact import binfile
from pynpact import capproc
from pynpact.util import Hasher, reducedict, mkstemp_overwrite


logger = logging.getLogger('pynpact.steps.extract')
statuslog = logging.getLogger('pynpact.statuslog')


KEYS = ['GeneDescriptorKey1', 'GeneDescriptorKey2',
        'GeneDescriptorSkip1', 'GeneDescriptorSkip2',
        'base_file']


class ExtractStep(BaseStep):
    def enqueue(self, config):
        """Go through the genbank record pulling out gene names and locations
        $ extract MYCGE.gbk 0 gene 0 locus_tag > MYCGE.genes
        """
        config, hash = get_hash(config)
        gbkfile = config['base_file']
        target_file = self.derive_filename(gbkfile, hash, 'genes')
        return self.executor.enqueue(
            delay(_extract)(gbkfile, target_file, config),
            id=target_file)


def get_hash(config):
    config = reducedict(config, KEYS)
    h = Hasher()
    h.hashdict(config)
    h.hashfiletime(binfile('extract'))
    h.hashfiletime(config['base_file'])
    return config, h.hexdigest()


def _extract(gbkfile, target_file, config):
    with mkstemp_overwrite(target_file) as out:
        statuslog.debug(
            "Extracting genes in %s.", os.path.basename(gbkfile))
        cmd = [binfile("extract"), gbkfile,
               config['GeneDescriptorSkip1'], config['GeneDescriptorKey1'],
               config['GeneDescriptorSkip2'], config['GeneDescriptorKey2']]
        capproc.capturedCall(
            cmd, check=True, logger=logger, stdout=out, stderr=sys.stderr)
    return target_file
