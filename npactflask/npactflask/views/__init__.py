import flask
import os
import os.path
from flask import url_for, flash, redirect, safe_join
from werkzeug.exceptions import NotFound
from npactflask import app

logger = app.logger


def library_root():
    abspath = getabspath("library", False)
    if not os.path.exists(abspath):
        logger.info("Creating library root at %s", abspath)
        os.makedirs(abspath)
    return abspath


def getrelpath(abspath):
    return os.path.relpath(abspath, app.config['UPLOADS'])


def is_clean_path(path):
    return (not path.startswith('/')) and os.path.normpath(path) == path


class ImmediateHttpResponseException (Exception):
    httpResponse = None

    def __init__(self, httpResponse):
        self.httpResponse = httpResponse


class RedirectException(Exception):
    url = None

    def __init__(self, url):
        self.url = url


def getabspath(relpath, raise_on_missing=True):
    abspath = safe_join(app.config['UPLOADS'], relpath)
    if raise_on_missing and not os.path.exists(abspath):
        raise NotFound(relpath)
    return abspath


@app.context_processor
def vSTATIC():
    def STATICV(filename):
        if app.config['DEBUG']:
            vnum = os.path.getmtime(os.path.join(app.static_folder, filename))
        else:
            vnum = app.config['VERSION']
        return (url_for('static', filename=filename, vnum=vnum))
    return dict(vSTATIC=STATICV)


def index():
    return flask.render_template('start.html')


def about():
    return flask.render_template('about.html')


def view_none():
    flash(
        "No genome source selected, please upload one, "
        "or go to the library and select one.")
    return redirect(url_for('index', filename=''))
