#!/usr/bin/env python
"""This script scans for apparently unused files in the upload path
and deletes them.

'unused files' are considered to be ones for which the atime is older
than 14 days (see ATIME_DEFAULT)
"""

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
    logger.info("Cleaning older than: %d @ %r", days, path)

    cmd=["find", path, "-atime", "+" + str(days), "-delete"]
    return capproc.capturedCall(cmd, logger=logger, stdin=False, stderr_level=logging.WARNING)


if __name__ == '__main__' :
    #make sure django settings is setup; used for the path source and logging config
    if 'DJANGO_SETTINGS_MODULE' not in os.environ:
        os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'
    from django.conf import settings


    parser = OptionParser("""usage: %prog [options]

    This script scans for apparently unused files in the upload path
    and the taskqueue path and deletes them.

    'unused files' are considered to be ones for which the atime is
    older than X days; set in django settings or via the --atime flag.
    """)

    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")

    parser.add_option("-a", "--atime", action="store", dest="atime",
                      default=settings.ATIME_DEFAULT,
                      help="argument to find's atime predicate for how many "
                      "days since it has been accessed before we decide to "
                      "delete it. Defaults to %s" % settings.ATIME_DEFAULT)

    (options,args) = parser.parse_args()


    if options.verbose:
        #logger is set to WARNING by default
        logger.setLevel(logging.DEBUG)

    try:
        days = int(options.atime)
    except:
        logger.error("Invalid number of days to keep for; must be an integer.")
        sys.exit(1)

    try:
        rc1 = clean(path(settings.MEDIA_ROOT).realpath(), days, options.verbose)
        rc2 = clean(path(settings.TQ_DIR).realpath(), days, options.verbose)
        sys.exit(rc1 or rc2)
    except SystemExit:
        raise
    except:
        logger.exception("Error during cleanup.")
        sys.exit(1)
