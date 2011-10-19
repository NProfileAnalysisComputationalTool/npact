#!/usr/bin/env python

ATIME_DEFAULT=14

import os
import logging
import subprocess

from optparse import OptionParser
from path import path

from pynpact import util
from spatweb import library_root


logger = logging.getLogger('cleanup')



def clean(path, days, verbose=False):
    logger.info("Cleaning up; days:%d, %r", days, path)

    if verbose:
        cmd=["find", path, "-atime", "+" + str(days),
             "-exec", "rm","-v", "{}", "+"]
    else:
        cmd=["find", path, "-atime", "+" + str(days),
             "-exec", "rm", "{}", "+"]
    util.capturedCall(cmd,logger=logger,stderr_level=logging.WARNING)


def setup_logger(verbose):
    logging.basicConfig(level=(verbose and logging.DEBUG or logging.WARNING),
                        format="%(asctime)s %(name)-10s %(levelname)-8s %(message)s",
                        datefmt='%H:%M:%S')
    if verbose:
        logging.getLogger('').setLevel(logging.DEBUG)

if __name__ == '__main__' :

    parser = OptionParser("""usage: %prog [options] [path]
        src <- The file to decrypt. If a directory, decrypt the most recent in there.""")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")

    parser.add_option("-a", "--atime", action="store", dest="atime",
                      default=ATIME_DEFAULT,
                      help="""argument to find's atime predicate for how many days since it has been accessed before we decide to delete it. Defaults to %s""" % ATIME_DEFAULT)

    (options,args) = parser.parse_args()
    setup_logger(options.verbose)
    
    clean(path(__file__).dirname().joinpath('../webroot/uploads').realpath(),
          options.atime,
          options.verbose)
