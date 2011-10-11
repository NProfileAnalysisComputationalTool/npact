#!/usr/bin/env python

import logging
import subprocess
from optparse import OptionParser


from django.conf import settings

from spat import library_root


logger = logging.getLogger(__name__)



def clean(path, days):
    logger.info("Cleaning up; days:%d, path:%r", days, path)
    
    cmd=["find", path, "-atime", "+" + days,
         "-exec", "rm", "{}", "+"]


def setup_logger(verbose):
    logging.basicConfig(level=(verbose and logging.DEBUG or logging.INFO),
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
                      default=settings.settings.MEDIA_RETAIN_FOR,
                      help="""argument to find's atime predicate for how many days since it has been accessed before we decide to delete it.""")

    (options,args) = parser.parse_args()
    setup_logger(options.verbose)

    
    #parser.add_option("-n", "--dry-run", action)
