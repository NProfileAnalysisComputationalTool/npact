from npactflask import app
import os
import logging
from logging.handlers import WatchedFileHandler

logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)
fmt = "%(asctime)s %(levelname)-7s %(name)-20s| %(message)s"
logdir = app.config['WEBROOT'] / 'logs'
logdir.makedirs_p()
fh = WatchedFileHandler(logdir / 'main.log')
fh.setFormatter(logging.Formatter(fmt))
logger.addHandler(fh)

if os.environ.get('ENV', 'dev').lower() in ('dev', 'development'):
    # Causes duplicate logging
    # app.config['DEBUG'] = True
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter(fmt, datefmt='%H:%M:%S'))
    logger.addHandler(sh)
    app.logger.handlers = []
