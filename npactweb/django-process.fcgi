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
import os, os.path, sys

os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'  # or whatever

from flup.server.fcgi import WSGIServer
from django.core.handlers.wsgi import WSGIHandler


WSGIServer(WSGIHandler(), debug=True).run()
