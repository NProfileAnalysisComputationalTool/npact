#!/usr/bin/env python

#fcgi entry point to the application


#this is supposed to work, if it doesn't remember this address:
# http://stackoverflow.com/questions/527237/unhandled-exception-in-flup
import os

os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'  # or whatever

from flup.server.fcgi import WSGIServer
from django.core.handlers.wsgi import WSGIHandler
from django.conf import settings

WSGIServer(WSGIHandler(), debug=settings.DEBUG).run()
