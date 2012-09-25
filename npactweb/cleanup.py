#!/usr/bin/env python
"""This script scans for apparently unused files in the upload path
and deletes them.

'unused files' are considered to be ones for which the atime is older
than 14 days (see ATIME_DEFAULT)
"""

#Default number of days that files can remain unaccessed before deletion.
#nb: this setting is duplicated in settings module for printing on the website.
ATIME_DEFAULT=14

import os
import logging
import subprocess
import sys
from optparse import OptionParser

from path import path

from pynpact import capproc
from npactweb import library_root


logger = logging.getLogger('cleanup')



def clean(path, days, verbose=False):
    "This function actually runs the find and delete."
    logger.info("Cleaning up; days:%d, %r", days, path)

    if verbose:
        cmd=["find", path, "-atime", "+" + str(days),
             "-exec", "rm", "-f", "-v", "{}", "+"]
    else:
        cmd=["find", path, "-atime", "+" + str(days),
             "-exec", "rm", "-f", "{}", "+"]
    return capproc.capturedCall(cmd, logger=logger, stdin=False, stderr_level=logging.WARNING)


def setup_logger(verbose):
    logging.basicConfig(level=(verbose and logging.DEBUG or logging.WARNING),
                        format="%(asctime)s %(name)-10s %(levelname)-8s %(message)s",
                        datefmt='%H:%M:%S')
    if verbose:
        logging.getLogger('').setLevel(logging.DEBUG)

if __name__ == '__main__' :

    parser = OptionParser("""usage: %prog [options]

    This script scans for apparently unused files in the upload path
    and deletes them.

    'unused files' are considered to be ones for which the atime is
    older than 14 days (see ATIME_DEFAULT)
    """)

    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")

    parser.add_option("-a", "--atime", action="store", dest="atime",
                      default=ATIME_DEFAULT,
                      help="""argument to find's atime predicate for how many days since it has been accessed before we decide to delete it. Defaults to %s""" % ATIME_DEFAULT)

    (options,args) = parser.parse_args()
    setup_logger(options.verbose)

    if not 'DJANGO_SETTINGS_MODULE' in os.environ:
        os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

    try:
        from django.conf import settings
        path = path(settings.MEDIA_ROOT).realpath()
        rc = clean(path, int(options.atime), options.verbose)
        sys.exit(rc)
    except SystemExit:
        raise
    except:
        logger.exception("Error during cleanup.")
        sys.exit(1)
