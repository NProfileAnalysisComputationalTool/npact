from gevent.monkey import patch_all
patch_all(subprocess=True)

import pkg_resources
from flask import Flask, redirect
from werkzeug.wsgi import DispatcherMiddleware
from flask_mail import Mail

app = Flask(__name__)
app.config.from_object('npactflask.settings')
mail = Mail(app)

app.config['APPLICATION_ROOT'] = '/npact'
app.config['VERSION'] = pkg_resources.require('npactflask')[0].version


from npactflask import helpers, setuplogging
from npactflask.views import raw, about
from npactflask.views import run, start, management

SILENCE_UNUSED_WARNING = (helpers, setuplogging)


if app.config.get('APPLICATION_ROOT', '/') != '/':
    # ***  Handle running at /npact/  ***
    application = Flask('redirectapp')

    @application.route('/')
    def doredirect():
        return redirect(app.config['APPLICATION_ROOT'])

    # Load a redirect app at the root URL to redirect us to the target app.
    # Serve app at APPLICATION_ROOT for localhost development.
    application = DispatcherMiddleware(application, {
        app.config['APPLICATION_ROOT']: app,
    })
else:
    application = app
