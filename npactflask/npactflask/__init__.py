import pkg_resources
from flask import Flask, redirect
from werkzeug.wsgi import DispatcherMiddleware


app = Flask(__name__)
app.config['APPLICATION_ROOT'] = '/npact'
app.config['VERSION'] = pkg_resources.require('npactflask')[0].version
app.secret_key = 'cf0cb53d-1ff1-4074-a65d-977831de66af'

from npactflask import settings, helpers
from npactflask.views import raw, about
from npactflask.views import run, start, management

app.add_url_rule('/', 'start', view_func=start.view, methods=['POST', 'GET'])

app.add_url_rule('/run/<path:path>', 'run_frame', view_func=run.run_frame)
app.add_url_rule('/about', view_func=about)
app.add_url_rule('/runstatus/<path:path>', view_func=run.run_status)
app.add_url_rule('/kickstart/<path:path>', view_func=run.kickstart)
app.add_url_rule('/translate', view_func=run.translate, methods=['POST'])
app.add_url_rule('/acgt_gamma_file_list/<path:path>',
                 view_func=run.acgt_gamma_file_list)

app.add_url_rule('/efetch/<int:id>', view_func=start.efetch)
app.add_url_rule('/management', 'management', view_func=management.view,
                 methods=['POST', 'GET'])

app.add_url_rule('/raw/<path:path>', 'raw', raw)



redirectapp = Flask('redirectapp')


@redirectapp.route('/')
def doredirect():
    return redirect(app.config['APPLICATION_ROOT'])

# Load a redirect app at the root URL to redirect us to the target app.
# Serve app at APPLICATION_ROOT for localhost development.
app_with_redirect = DispatcherMiddleware(redirectapp, {
    app.config['APPLICATION_ROOT']: app,
})
