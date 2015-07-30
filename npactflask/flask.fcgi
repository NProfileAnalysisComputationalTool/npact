#!/usr/bin/env python
#http://flask.pocoo.org/docs/0.10/deploying/fastcgi/#creating-a-fcgi-file
from logging.handlers import WatchedFileHandler
import logging
from logging import getLogger, Formatter

from flask import Flask, redirect

from flup.server.fcgi import WSGIServer
from npactflask import app, app_with_redirect
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

application = app_with_redirect
if __name__ == '__main__':
    WSGIServer(application).run()
