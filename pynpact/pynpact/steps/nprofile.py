"""Do the CG ratio calculations.

        $ CG MYCGE.gbk 1 580074 201 51 3 > MYCGE.CG200
        """
from __future__ import absolute_import
import logging
import sys
import json

from path import path
from pynpact import binfile
from pynpact import capproc, parsing
from pynpact.util import Hasher, reducedict
from pynpact.steps import producer, enqueue


logger = logging.getLogger('pynpact.steps.nprofile')
statuslog = logging.getLogger('pynpact.statuslog')


BIN = binfile('nprofile')

KEYS = ['nucleotides', 'length', 'window_size', 'step', 'period', 'ddna']
OUTPUTKEY = 'File_list_of_nucleotides_in_200bp windows'
JSONOUTPUTKEY = 'nprofileData'


def plan(config, executor):
    if 'nprofile' in config:
        return
    config['nprofile'] = True

    parsing.length(config)
    rconfig = reducedict(config, KEYS)
    h = Hasher()
    h.hashdict(rconfig)
    h.hashfiletime(BIN)
    hash = h.hexdigest()
    target = parsing.derive_filename(config, hash, 'nprofile')
    config[OUTPUTKEY] = target
    config[JSONOUTPUTKEY] = target + '.json'
    jobs = enqueue(_nprofile, executor, rconfig, target)
    enqueue(_nprofile_to_json, executor, {OUTPUTKEY: target},
            config[JSONOUTPUTKEY], after=jobs)
    return jobs


@producer()
def _nprofile(config, out):
    statuslog.info("Calculating n-profile.")
    cmd = [BIN, '-b', ''.join(config["nucleotides"]).upper(),
           config['ddna'], 1, config['length'],
           config['window_size'], config['step'], config['period']]
    capproc.capturedCall(
        cmd, stdout=out, stderr=sys.stderr,
        logger=logger, check=True)


def parse_nprofile(ifile):
    # Read the file in and split the fields on space
    nprofile_lines = [l.split() for l in path(ifile).lines(retain=False)]
    # convert the strings to numbers
    data = [(int(c), float(x), float(y), float(z))
            for (c, x, y, z) in nprofile_lines]
    return data


@producer()
def _nprofile_to_json(config, out):
    keys = ['coordinate', 'r', 'g', 'b']
    data = parse_nprofile(config[OUTPUTKEY])
    data = [dict(zip(keys, d)) for d in data]
    json.dump(data, out)
