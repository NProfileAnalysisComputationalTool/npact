#!/usr/bin/env python
"""This script scans for apparently unused files in the upload path
and deletes them.

'unused files' are considered to be ones for which the atime is older
than 14 days (see ATIME_DEFAULT)
"""
import logging
import sys
from optparse import OptionParser


logger = logging.getLogger('cleanup')


if __name__ == '__main__':
    from npactflask import app, cleanup

    parser = OptionParser("""usage: %prog [options]

    This script scans for apparently unused files in the upload path
    and the taskqueue path and deletes them.

    'unused files' are considered to be ones for which the atime is
    older than X days; set in django settings or via the --atime flag.
    """)

    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")

    parser.add_option("-a", "--atime", action="store", dest="atime",
                      default=app.config['ATIME_DEFAULT'],
                      help="argument to find's atime predicate for how many "
                      "days since it has been accessed before we decide to "
                      "delete it. Defaults to %default")

    (options, args) = parser.parse_args()
    if options.verbose:
        logger.setLevel(logging.DEBUG)

    try:
        days = int(options.atime)
    except:
        logger.error("Invalid number of days to keep for; must be an integer.")
        sys.exit(1)

    try:
        if cleanup.cleanup_old_files(days):
            logger.info("Success!")
    except SystemExit:
        raise
    except:
        logger.exception("Error during cleanup.")
        sys.exit(1)
