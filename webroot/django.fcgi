#!../ve/bin/python

#this is supposed to work, if it doesn't remember this address:
# http://stackoverflow.com/questions/527237/unhandled-exception-in-flup
import os, os.path, sys

#move up a directory so that the 'import settings' will work. We don't
#just put this file in that directory b/c we want to limit *Apache* to
#the webroot directory
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))


from flup.server.fcgi import WSGIServer

os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'  # or whatever

from django.core.handlers.wsgi import WSGIHandler


WSGIServer(WSGIHandler(), debug=True).run()
