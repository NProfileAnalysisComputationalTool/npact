#!/usr/bin/env python
"""

# This should be a daemon or a cron task?

# there will be a queue directory that holds to be done tasks. 

#Do tasks in the directory in date created order.

# ENQUEUE:  mkstemp should guarantee uniquenames

# DEQUEUE: list files in directory, sort by date modified, start working on job.
#  * multiple pool processes?
#   * have one process doing the dequeue, use multiprocessing module beyond that.
 


"""

import os
import logging
import tempfile
import sys
from optparse import OptionParser

from path import path

from pynpact import capproc
from npactweb import library_root

if not 'DJANGO_SETTINGS_MODULE' in os.environ:
    os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'
from django.conf import settings


logger = logging.getLogger('taskqueue')

def get_queue_directory():
    return path(settings.QUEUE_DIR).realpath()

def enqueue(task, args):
    """Enqueue a task. 

    Task should be the name of afully referenced python callable and
    will be passed the given args. Everything should be picklable.
    """
    fd,abspath = tempfile.mkstemp(dir=get_queue_directory())
    #write task as first line
    #write args as pickled rest
    

def status(id):
    """Gets the current status of the identified task"""
    pass


def dequeue():
    pass

def setup_logger(verbose):
    logging.basicConfig("", level=(verbose and logging.DEBUG or logging.WARNING),
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

    if not 'DJANGO_SETTINGS_MODULE' in os.environ:
        os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

    try:
        from django.conf import settings
        
        sys.exit(rc)
    except SystemExit:
        raise
    except:
        logger.exception("Error during cleanup.")
        sys.exit(1)
