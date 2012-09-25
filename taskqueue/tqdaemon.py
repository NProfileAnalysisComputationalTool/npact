#!/usr/bin/env python
"""

# This should be a daemon or a cron task?  daemon i think


http://pypi.python.org/pypi/python-daemon/
"""
import logging
import multiprocessing
import sys
from optparse import OptionParser


logger = logging.getLogger('taskqueue.daemon')

def setup_logger(verbose):
    tq_logger = logging.getLogger('')
    tq_logger.propagate=False
    tq_logger.setLevel(logging.DEBUG if verbose else logging.WARNING)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter("%(asctime)s %(processName)s/%(name)s "
                                           "%(levelname)-8s %(message)s",
                                           datefmt='%H:%M:%S'))
    tq_logger.addHandler(handler)
    mp_logger = multiprocessing.get_logger()
    mp_logger.addHandler(handler)
    mp_logger.setLevel(logging.INFO)



if __name__ == '__main__' :
    parser = OptionParser("""usage: %prog [options]

    """)

    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")

    (options,args) = parser.parse_args()
    setup_logger(options.verbose)
    try:
        from taskqueue import server
        server.start_everything()
        sys.exit(0)
    except SystemExit:
        raise
    except:
        logger.exception("Error in daemon")
        sys.exit(1)
