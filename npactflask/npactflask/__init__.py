import pkg_resources
from flask import Flask, send_from_directory


app = Flask(__name__)
app.config['APPLICATION_ROOT'] = '/npact'
app.config['VERSION'] = pkg_resources.require('npactflask')[0].version
app.secret_key = 'cf0cb53d-1ff1-4074-a65d-977831de66af'


from npactflask import settings
from npactflask import views
from npactflask.views import run, start, management


app.add_url_rule('/', view_func=start.view, methods=['POST', 'GET'],
                 endpoint="start")
app.add_url_rule('/run/<path:path>', view_func=run.run_frame)
app.add_url_rule('/about', view_func=views.about)
app.add_url_rule('/runstatus/<path:path>', view_func=run.run_status)
app.add_url_rule('/kickstart/<path:path>', view_func=run.kickstart)
app.add_url_rule('/translate', view_func=run.translate, methods=['POST'])
app.add_url_rule('/acgt_gamma_file_list/<path:path>',
                 view_func=run.acgt_gamma_file_list)
app.add_url_rule('/efetch/<int:id>', view_func=start.efetch)
app.add_url_rule(r'/^(run|config)^', view_func=views.view_none)
app.add_url_rule('/management', view_func=management.view,
                 methods=['POST', 'GET'], endpoint='management')
app.add_url_rule('/file-upload', view_func=start.file_upload, methods="POST")


@app.route('/raw/<path:path>')
def raw(path):
    return send_from_directory(app.config['UPLOADS'], path)
