import os.path
import sys
import logging
import errno
import tempfile
import shutil
from contextlib import contextmanager
from pynpact import capproc, parsing
from pynpact import binfile, DATAPATH
from pynpact.util import \
    Hasher, reducedict, delay
log = logging.getLogger(__name__)

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

    executor.enqueue(delay(_acgt_gamma)(rconfig, outdir),
                     tid=outdir)

    log.debug("Adding prediction filenames to config dict.")
    # strip 4 characters off here b/c that's how acgt_gamma does
    # it at about lines 262-270
    j = lambda ext: os.path.join(
        outdir, os.path.basename(config['filename'])[:-4] + ext)
    config['File_of_new_CDSs'] = j(".newcds")
    config['File_of_published_rejected_CDSs'] = j(".modified")
    config['File_of_G+C_coding_potential_regions'] = j('.profiles')
    config[OUTPUTKEY] = outdir
    return outdir


def _acgt_gamma(config, target_dir):
    if os.path.exists(target_dir):
        # all timestamps are hashed in at this point, if anything
        # changes it has a new name. Hence if this dir exists we're
        # already done.
        return

    cmd = [BIN, "-q", config['filename']]
    if 'significance' in config:
        cmd.append(config['significance'])

    # TODO: acgt actually takes the string of characters to skip,
    # not the length.
    with mkdtemp() as dtemp:
        #statuslog.info("Identifying ORFs with 3-base compositional periodicity with significance: %s.", config.get('significance'))
        log.debug("Starting prediction program in %s", dtemp)
        capproc.capturedCall(
            cmd, cwd=dtemp, check=True,
            env={'BASE_DIR_THRESHOLD_TABLES': DATAPATH},
            stderr=sys.stderr,
            logger=log)
        log.debug("Renaming from %s to %s", dtemp, target_dir)
        try:
            os.rename(dtemp, target_dir)
        except OSError as e:
            if e.errno == errno.EEXIST:
                shutil.rmtree(dtemp)
            else:
                raise
    return target_dir


@contextmanager
def mkdtemp(**kwargs):
    """A wrapper for tempfile.mkdtemp that always cleans up.

    This wrapper sets defaults based on the class values."""
    path = tempfile.mkdtemp(**kwargs)
    try:
        yield path
    finally:
        log.debug("Cleaning up mkdtemp %r", path)
        shutil.rmtree(path, ignore_errors=True)
