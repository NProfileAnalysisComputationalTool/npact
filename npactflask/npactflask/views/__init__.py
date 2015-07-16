import flask
import logging
import os
import os.path
from flask import url_for, send_from_directory, flash, redirect
from werkzeug.exceptions import NotFound
from npactflask import app
from npactflask import settings

logger = logging.getLogger(__name__)


def library_root():
    abspath = getabspath("library", False)
    if not os.path.exists(abspath):
        logger.info("Creating library root at %s", abspath)
        os.makedirs(abspath)
    return abspath


def getrelpath(abspath):
    return os.path.relpath(abspath, settings.MEDIA_ROOT)


def is_clean_path(path):
    return (not path.startswith('/')) and os.path.normpath(path) == path


def assert_clean_path(path, request,
                      message="Path contained illegal characters. Please "
                      "select a new GBK file.",
                      destination='start'):
    if not is_clean_path(path):
        flash(message)
        raise RedirectException(url_for(destination))


class ImmediateHttpResponseException (Exception):
    httpResponse = None

    def __init__(self, httpResponse):
        self.httpResponse = httpResponse


class RedirectException(Exception):
    url = None

    def __init__(self, url):
        self.url = url


def getabspath(relpath, raise_on_missing=True):
    if not is_clean_path(relpath):
        logger.error("Illegal path submitted", relpath)
        raise NotFound("Bad characters")
    abspath = os.path.join(settings.MEDIA_ROOT, relpath)
    if raise_on_missing and not os.path.exists(abspath):
        raise NotFound
    return abspath


@app.context_processor
def STATIC():
    return dict(STATIC_URL=url_for('static', filename=''))


def index():
    return flask.render_template('start.html')


def about():
    return flask.render_template('about.html')


def view_none():
    flash(
        "No genome source selected, please upload one, "
        "or go to the library and select one.")
    return redirect(url_for('index', filename=''))


def static_serve_wrapper(path):
    if app.debug:
        return send_from_directory(settings.MEDIA_ROOT, path)
