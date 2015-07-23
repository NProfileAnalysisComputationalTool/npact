#!/usr/bin/env python
#http://flask.pocoo.org/docs/0.10/deploying/fastcgi/#creating-a-fcgi-file
from logging.handlers import WatchedFileHandler
import logging
from logging import getLogger, Formatter, StreamHandler

from flask import Flask, redirect
from werkzeug.wsgi import DispatcherMiddleware
from flup.server.fcgi import WSGIServer
from npactflask import app
from npactflask.settings import ppath
from taskqueue.tqdaemon import tqdaemonlog



app.config['DEBUG'] = False
tqdaemonlog()
vFormatter = Formatter('%(asctime)s %(process)d:%(thread)d %(name)-15s'
                       ' %(levelname)-7s| %(message)s')
root = getLogger('')
roothandler = WatchedFileHandler(ppath('logs') / 'main.log')
roothandler.setFormatter(vFormatter)
root.addHandler(roothandler)
root.setLevel(logging.INFO)

redirectapp = Flask('redirectapp')


@redirectapp.route('/')
def doredirect():
    return redirect(app.config['APPLICATION_ROOT'])

# Load a redirect app at the root URL to redirect us to the target app.
# Serve app at APPLICATION_ROOT for localhost development.
application = DispatcherMiddleware(redirectapp, {
    app.config['APPLICATION_ROOT']: app,
})

if __name__ == '__main__':

    #TODO: this probably needs the DispatcherMiddleware similar to runserver
    WSGIServer(application).run()
