#!/usr/bin/env python
"""

# This should be a daemon or a cron task?  daemon i think


http://pypi.python.org/pypi/python-daemon/
"""
import logging
from optparse import OptionParser

logger = logging.getLogger('taskqueue.tqdaemon')


if __name__ == '__main__':
    parser = OptionParser("""usage: %prog [options]

    """)

    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")

    (options, args) = parser.parse_args()
    import taskqueue.tqdaemon
    taskqueue.setup_logger(options.verbose)
    taskqueue.tqdaemon.daemonize()
