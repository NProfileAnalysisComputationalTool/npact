#!/usr/bin/env python

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
