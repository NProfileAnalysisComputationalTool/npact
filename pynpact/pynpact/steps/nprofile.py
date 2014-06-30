"""Do the CG ratio calculations.

        $ CG MYCGE.gbk 1 580074 201 51 3 > MYCGE.CG200
        """
from __future__ import absolute_import
import logging
import os.path
import sys

from pynpact.steps import BaseStep, delay
from pynpact import binfile
from pynpact import capproc
from pynpact.util import Hasher, reducedict, mkstemp_overwrite


logger = logging.getLogger('pynpact.steps.nprofile')
statuslog = logging.getLogger('pynpact.statuslog')


BIN = binfile('nprofile')

KEYS = ['nucleotides', 'length', 'window_size', 'step', 'period', 'base_file']


class NprofileStep(BaseStep):
    def enqueue(self):
        config = reducedict(self.config, KEYS)
        base_file = config['base_file']
        h = Hasher()
        h.hashdict(config)
        h.hashfiletime(BIN)
        h.hashfiletime(base_file)
        hash = h.hexdigest()
        filename = self.derive_filename(base_file, hash, 'nprofile')
        self.executor.enqueue(
            delay(_nprofile)(base_file, filename, config),
            id=filename)
        return filename


def _nprofile(base_file, target_file, config):
    statuslog.info("Calculating n-profile.")
    cmd = [
        BIN, '-b', ''.join(config["nucleotides"]),
        base_file, 1, config['length'],
        config['window_size'], config['step'], config['period']]
    with mkstemp_overwrite(target_file) as out:
        capproc.capturedCall(
            cmd, stdout=out, stderr=sys.stderr,
            logger=logger, check=True)
