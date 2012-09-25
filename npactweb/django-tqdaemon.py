#!/usr/bin/env python

# We have two copies of this file to correspond with two process pools
# managed by apache. Apache requires a different target in order to
# start a second process pool (class).
#
# There are two process pools:
#
# * One handles general website traffic that is expected to complete
#   quickly
# * One handles all calls to the underlying C which could potentially
#   be slow to respond.


#this is supposed to work, if it doesn't remember this address:
# http://stackoverflow.com/questions/527237/unhandled-exception-in-flup
import os
import sys
import logging
import taskqueue

os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'  # or whatever


if __name__ == '__main__' :
    from optparse import OptionParser
    parser = OptionParser("""usage: %prog [options]

    """)

    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")

    (options,args) = parser.parse_args()
    
    try:
        #import and read from settings to get it to configure logging
        from django.conf import settings
        settings.LOGGING 

            
    except Exception, e:
        sys.stderr.write('error configuring logging: ' + str(e))

    taskqueue.daemonize()
