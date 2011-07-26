import logging, os.path, tempfile, time, shutil
from optparse import OptionParser
from contextlib import contextmanager

from Bio import SeqIO

import util

logger = logging.getLogger(__name__)

def reduce_genbank(gbpath) :
    def filterfun(outfile) :
        with open(gbpath,'r') as infile:
            for l in infile :
                outfile.write(l)
                if l.startswith("ORIGIN") :
                    outfile.write("//\n")
                    return

    return util.safe_produce_new(gbpath,".noseq", filterfun,
                                 replace_ext=False, logger=logger)
