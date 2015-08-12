#!/usr/bin/env python

from gevent.monkey import patch_all
patch_all(subprocess=True)

import logging
logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)
fmt = "%(asctime)s %(levelname)-7s %(name)-20s| %(message)s"
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter(fmt, datefmt='%H:%M:%S'))
logger.addHandler(sh)

from npactflask import app, app_with_redirect

app.config['DEBUG'] = True



application = app_with_redirect


if __name__ == '__main__':
    # Relevant documents:
    # http://werkzeug.pocoo.org/docs/middlewares/
    # http://flask.pocoo.org/docs/patterns/appdispatch/
    from werkzeug.serving import run_simple
    run_simple('127.0.0.1', 5000, application,
               use_reloader=True, use_debugger=True, use_evalex=True)
