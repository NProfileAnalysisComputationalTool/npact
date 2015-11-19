
import pkg_resources
import os
import logging
from logging.handlers import WatchedFileHandler
from flask import Flask, redirect
from flask_mail import Mail

app = Flask(__name__)
app.config.from_object('npactflask.settings')
mail = Mail(app)

app.config['APPLICATION_ROOT'] = '/npact'
app.config['VERSION'] = pkg_resources.require('npactflask')[0].version


logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)
fmt = "%(asctime)s %(levelname)-7s %(name)-20s| %(message)s"
logdir = app.config['LOGDIR'].makedirs_p()
fh = WatchedFileHandler(logdir / 'main.log')
fh.setFormatter(logging.Formatter(fmt))
logger.addHandler(fh)

if app.config['DEBUG']:
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter(fmt, datefmt='%H:%M:%S'))
    logger.addHandler(sh)
    app.logger.handlers = []  # we've superceded flask's handler, get rid of it


# all of the following register stuff in the environment
from npactflask import helpers
from npactflask.views import raw, about
from npactflask.views import run, start, management

application = app

 # 'HTTP_X_FORWARDED_FOR': '172.17.42.1',
 # 'HTTP_X_FORWARDED_HOST': '172.17.0.6',
 # 'HTTP_X_FORWARDED_SERVER': '172.17.0.6',
 # 'PATH_INFO': '/run/library/NC_007760.gbk',
 # 'QUERY_STRING': 'nucleotides=C&nucleotides=G&significance=0.01&startBase=0&endBase=5013479&basesPerGraph=15000&offset=0&mycoplasma=false',
 # 'REMOTE_ADDR': '127.0.0.1',
