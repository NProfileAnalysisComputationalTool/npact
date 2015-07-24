import os.path

from flask import url_for

from npactflask import app


# TODO: I think this is more simply a template_global:
#  http://flask.pocoo.org/docs/0.10/api/#flask.Flask.template_global
@app.context_processor
def vSTATIC():
    def STATICV(filename):
        if app.config['DEBUG']:
            vnum = os.path.getmtime(os.path.join(app.static_folder, filename))
        else:
            vnum = app.config['VERSION']
        return (url_for('static', filename=filename, vnum=vnum))
    return dict(vSTATIC=STATICV)
