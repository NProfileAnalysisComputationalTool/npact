#!/usr/bin/env python

from gevent.monkey import patch_all
patch_all(subprocess=True)

from logging.handlers import WatchedFileHandler
import logging
from logging import getLogger, Formatter


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
