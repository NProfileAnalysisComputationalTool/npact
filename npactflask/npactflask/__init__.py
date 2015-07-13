from flask import Flask

app = Flask(__name__)
app.config['APPLICATION_ROOT'] = '/npact'


from npactflask import views


app.add_url_rule('/', view_func=views.index)
app.add_url_rule('/run/<path:path>', view_func=views.run_frame)
app.add_url_rule('/about', view_func=views.about)
app.add_url_rule('/runstatus/<path:path>', view_func=views.run_status)
app.add_url_rule('/kickstart/<path:path>', view_func=views.kickstart)
app.add_url_rule('/translate', view_func=views.translate)
app.add_url_rule('/acgt_gamma_file_list/<path:path>', view_func=views.acgt_gamma_file_list)
app.add_url_rule('/efetch/<int:path>', view_func=views.efetch)
app.add_url_rule(r'/^(run|config)^', view_func=views.view_none)
app.add_url_rule('/raw/<path:path>', view_func=views.static_serve_wrapper, endpoint='raw')
app.add_url_rule('/management', view_func=views.index)
