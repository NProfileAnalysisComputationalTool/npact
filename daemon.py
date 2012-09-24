#!/usr/bin/env python
"""

# This should be a daemon or a cron task?  daemon i think


http://pypi.python.org/pypi/python-daemon/
"""

import logging
import sys
from optparse import OptionParser

from path import path


logger = logging.getLogger('taskqueue.daemon')


def setup_logger(verbose):
    logging.basicConfig( level=(verbose and logging.DEBUG or logging.WARNING),
                        format="%(asctime)s %(name)-10s %(levelname)-8s %(message)s",
                        datefmt='%H:%M:%S')
    if verbose:
        logging.getLogger('').setLevel(logging.DEBUG)


if __name__ == '__main__' :
    parser = OptionParser("""usage: %prog [options]
    
    """)

    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")

    (options,args) = parser.parse_args()
    setup_logger(options.verbose)
    try:
        from taskqueue import server
        server.Server().run()
        sys.exit(0)
    except SystemExit:
        raise
    except:
        logger.exception("Error during cleanup.")
        sys.exit(1)
