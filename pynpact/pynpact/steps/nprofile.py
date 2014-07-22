"""Do the CG ratio calculations.

        $ CG MYCGE.gbk 1 580074 201 51 3 > MYCGE.CG200
        """
from __future__ import absolute_import
import logging
import sys


from pynpact import binfile
from pynpact import capproc, parsing
from pynpact.util import Hasher, reducedict, mkstemp_overwrite, delay


logger = logging.getLogger('pynpact.steps.nprofile')
statuslog = logging.getLogger('pynpact.statuslog')


BIN = binfile('nprofile')

KEYS = ['nucleotides', 'length', 'window_size', 'step', 'period', 'filename']
OUTPUTKEY = 'File_list_of_nucleotides_in_200bp windows'


def plan(config, executor):
    if 'nprofile' in config:
        return
    config['nprofile'] = True

    parsing.length(config)
    rconfig = reducedict(config, KEYS)
    h = Hasher()
    h.hashdict(rconfig)
    h.hashfiletime(BIN)
    h.hashfiletime(config['filename'])
    hash = h.hexdigest()
    target = parsing.derive_filename(config, hash, 'nprofile')
    executor.enqueue(
        delay(_nprofile)(rconfig, target),
        tid=target)
    config[OUTPUTKEY] = target
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
