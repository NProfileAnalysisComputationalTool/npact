#!/usr/bin/env python
"""

# This should be a daemon or a cron task?  daemon i think


http://pypi.python.org/pypi/python-daemon/
"""

import logging
import sys
from optparse import OptionParser

import taskqueue



logger = logging.getLogger('taskqueue.daemon')



if __name__ == '__main__' :
    parser = OptionParser("""usage: %prog [options]

    """)

    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")

    (options,args) = parser.parse_args()
    taskqueue.setup_logger(options.verbose)
    try:
        from taskqueue import server
        server.start_everything()
        sys.exit(0)
    except SystemExit:
        raise
    except:
        logger.exception("Error in daemon")
        sys.exit(1)
