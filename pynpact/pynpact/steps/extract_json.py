"""Go through the genbank record pulling out CDS as json"""
from __future__ import absolute_import
import logging
import os.path
from pynpact import genbank, parsing
from pynpact.util import Hasher, reducedict
from pynpact.steps import producer, enqueue


logger = logging.getLogger('pynpact.steps.extract_json')
statuslog = logging.getLogger('pynpact.statuslog')
KEYS = ['filename', 'stderr']


def plan(config, executor):
    if parsing.isgbk(config):
        logger.debug(
            "GBK file, extracting known CDS to json %s", config['filename'])
        rconfig, hash = get_hash(config)
        target_file = parsing.derive_filename(config, hash, 'track.genes.json')
        config['InputCDSFileJson'] = target_file
        return enqueue(_extract, executor, rconfig, target_file)


def get_hash(config):
    config = reducedict(config, KEYS)
    h = Hasher()
    h.hashdict(config)
    h.hashfiletime(config['filename'])
    return config, h.hexdigest()


@producer()
def _extract(config, out):
    gbkfile = config['filename']
    statuslog.debug(
        "Extracting CDS to JSON %s.", os.path.basename(gbkfile))
    genbank.gbk_to_track_json(gbkfile, out)
    logger.debug("Finished extract_json")
