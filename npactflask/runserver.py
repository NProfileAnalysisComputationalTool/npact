#!/usr/bin/env python
import logging
from flask import Flask, redirect
from npactflask import app

logger = logging.getLogger('npactflask')
logger.setLevel(logging.DEBUG)
# logger.addHandler(logging.StreamHandler())
fmt = "%(asctime)s %(levelname)-7s %(name)-20s| %(message)s"
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter(fmt, datefmt='%H:%M:%S'))
logger.addHandler(sh)


redirectapp = Flask('redirectapp')


@redirectapp.route('/')
def doredirect():
    return redirect(app.config['APPLICATION_ROOT'])


if __name__ == '__main__':
    # Relevant documents:
    # http://werkzeug.pocoo.org/docs/middlewares/
    # http://flask.pocoo.org/docs/patterns/appdispatch/
    from werkzeug.serving import run_simple
    from werkzeug.wsgi import DispatcherMiddleware
    app.config['DEBUG'] = True
    # Load a redirect app at the root URL to redirect us to the target app.
    # Serve app at APPLICATION_ROOT for localhost development.
    application = DispatcherMiddleware(redirectapp, {
        app.config['APPLICATION_ROOT']: app,
    })
    run_simple('localhost', 5000, application,
               use_reloader=True, use_debugger=True, use_evalex=True)
