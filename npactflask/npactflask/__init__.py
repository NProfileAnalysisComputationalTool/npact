from gevent.monkey import patch_all
patch_all(subprocess=True)

import pkg_resources
from flask import Flask, redirect
from werkzeug.wsgi import DispatcherMiddleware

app = Flask(__name__)
app.config.from_object('npactflask.settings')
app.config['APPLICATION_ROOT'] = '/npact'
app.config['VERSION'] = pkg_resources.require('npactflask')[0].version


from npactflask import helpers, setuplogging
from npactflask.views import raw, about
from npactflask.views import run, start, management

SILENCE_UNUSED_WARNING = (helpers, setuplogging)


# ***  Handle running at /npact/  ***
redirectapp = Flask('redirectapp')



@redirectapp.route('/')
def doredirect():
    return redirect(app.config['APPLICATION_ROOT'])

# Load a redirect app at the root URL to redirect us to the target app.
# Serve app at APPLICATION_ROOT for localhost development.
app_with_redirect = DispatcherMiddleware(redirectapp, {
    app.config['APPLICATION_ROOT']: app,
})
