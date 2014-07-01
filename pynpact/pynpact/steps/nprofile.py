"""Do the CG ratio calculations.

        $ CG MYCGE.gbk 1 580074 201 51 3 > MYCGE.CG200
        """
from __future__ import absolute_import
import logging
import sys

from pynpact.steps import BaseStep
from pynpact import binfile
from pynpact import capproc
from pynpact.util import Hasher, reducedict, mkstemp_overwrite, delay


logger = logging.getLogger('pynpact.steps.nprofile')
statuslog = logging.getLogger('pynpact.statuslog')


BIN = binfile('nprofile')

KEYS = ['nucleotides', 'length', 'window_size', 'step', 'period', 'filename']


class NprofileStep(BaseStep):
    def enqueue(self):
        config = reducedict(self.config, KEYS)
        h = Hasher()
        h.hashdict(config)
        h.hashfiletime(BIN)
        h.hashfiletime(config['filename'])
        hash = h.hexdigest()
        target = self.derive_filename(hash, 'nprofile')
        self.executor.enqueue(
            delay(_nprofile)(config, target),
            id=target)
        return target


def _nprofile(config, target_file):
    statuslog.info("Calculating n-profile.")
    filename = config['filename']
    cmd = [
        BIN, '-b', ''.join(config["nucleotides"]),
        filename, 1, config['length'],
        config['window_size'], config['step'], config['period']]
    with mkstemp_overwrite(target_file) as out:
        capproc.capturedCall(
            cmd, stdout=out, stderr=sys.stderr,
            logger=logger, check=True)
