
from flask import url_for
from npactflask import app


@app.template_global()
def vSTATIC(filename):
    if app.config['DEBUG']:
        return url_for('static', filename=filename)
    else:
        return url_for('static',
                       filename=filename, vnum=app.config['VERSION'])
