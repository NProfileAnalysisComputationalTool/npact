import flask
import os
import os.path
from flask import safe_join, send_from_directory
from werkzeug.exceptions import NotFound
from npactflask import app

logger = app.logger


def raw(path):
    return send_from_directory(app.config['UPLOADS'], path)


def index():
    return flask.render_template('start.html')


def about():
    return flask.render_template('about.html')


def library_root():
    abspath = getabspath("library", False)
    if not os.path.exists(abspath):
        logger.info("Creating library root at %s", abspath)
        os.makedirs(abspath)
    return abspath


def getrelpath(abspath):
    return os.path.relpath(abspath, app.config['UPLOADS'])


def getabspath(relpath, raise_on_missing=True):
    abspath = safe_join(app.config['UPLOADS'], relpath)
    if raise_on_missing and not os.path.exists(abspath):
        raise NotFound(relpath)
    return abspath
