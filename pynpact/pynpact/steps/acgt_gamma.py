import os.path
import sys
import logging
from pynpact import capproc, parsing
from pynpact import binfile, DATAPATH
from pynpact.util import Hasher, reducedict, mkdtemp_rename
from pynpact.steps import producer, enqueue


log = logging.getLogger(__name__)
statuslog = logging.getLogger('pynpact.statuslog')

BIN = binfile("acgt_gamma")
OUTPUTKEY = 'acgt_gamma_output'


def plan(config, executor):
    "Identifying ORFs with significant 3-base periodicities."
    if config.get('skip_prediction', False):
        return

    assert os.path.exists(DATAPATH), \
        "Missing pynpact/data for acgt_gamma prediction. " \
        "Expected at " + DATAPATH

    rconfig = reducedict(config,
                         ['filename', 'significance', 'GeneDescriptorSkip1'])
    h = Hasher().hashdict(rconfig)
    h.hashfiletime(config['filename'])
    h.hashfiletime(BIN)
    outdir = parsing.derive_filename(config, h.hexdigest(), '.predict')

    log.debug("Adding prediction filenames to config dict.")
    # strip 4 characters off here b/c that's how acgt_gamma does
    # it at about lines 262-270
    j = lambda ext: os.path.join(
        outdir, os.path.basename(config['filename'])[:-4] + ext)
    config['File_of_new_CDSs'] = j(".newcds")
    config['File_of_published_rejected_CDSs'] = j(".modified")
    config['File_of_G+C_coding_potential_regions'] = j('.profiles')
    config[OUTPUTKEY] = outdir

    return enqueue(_acgt_gamma, executor, rconfig, outdir)


@producer(tmpmanager=mkdtemp_rename)
def _acgt_gamma(config, dtemp):
    cmd = [BIN, "-q", config['filename']]
    if 'significance' in config:
        cmd.append(config['significance'])

    statuslog.info(
        "Identifying ORFs with 3-base compositional periodicity "
        "with significance: %s.", config.get('significance'))
    log.debug("Starting prediction program in %s", dtemp)
    # TODO: acgt actually takes the string of characters to skip,
    # not the length.
    capproc.capturedCall(
        cmd, cwd=dtemp, check=True,
        env={'BASE_DIR_THRESHOLD_TABLES': DATAPATH},
        stderr=sys.stderr,
        logger=log)
