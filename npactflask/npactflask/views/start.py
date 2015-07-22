import flask
import logging
import os.path
import tempfile
import urllib2
from flask_wtf import Form
from wtforms import fields
from npactflask.views import settings, getrelpath
from npactflask.views import is_clean_path, library_root
from pynpact import util, entrez, parsing
from flask import url_for, redirect, request, flash, send_from_directory


logger = logging.getLogger(__name__)


def mksavefile(prefix):
    """Wrapper around creating the file to save uploaded files in.

    Returns the (fd, abspath, relpath)
"""
    # we're using tempfile to ensure the file we get is unique and
    # aren't overwriting another.
    fd, abspath = tempfile.mkstemp(
        dir=settings.MEDIA_ROOT, prefix=prefix)
    relpath = getrelpath(abspath)
    return (fd, abspath, relpath)


def get_ti(size):
    return fields.TextField(attrs={'size': size})


def file_upload():
    fileup = request.files.get('file')

    if not fileup:
        flash('Must upload a File')
        return redirect(url_for('start') + '?active=upload')
    logger.info("Checking Uploaded file %s", fileup)
    if not is_clean_path(fileup.name):
        flash("Illegal filename")
        return redirect(url_for('start') + '?active=upload')
    fd, savepath, relpath = mksavefile("up-%s-" % fileup.filename)
    logger.info("Saving uploaded file to %r", savepath)
    fileup.save(savepath)
    return redirect(url_for('run_frame', path=relpath))


def fetchurl():
    url = request.form.get('fetchurl')
    if not url:
        flash('URL Required in URL Field')
        return redirect(url_for('start', active='url'))
    logger.debug("Going to pull from %r", url)
    pull_req = urllib2.Request(url)
    if pull_req.get_type == "file":
        flash("Illegal URL")
        return redirect(url_for('start') + '?active=url')
    try:
        fh = urllib2.urlopen(pull_req)
        fd, savepath, relpath = mksavefile("url")
        util.stream_to_file(fh, savepath)
        return redirect(url_for('run_frame', path=relpath))
    except:
        logger.exception("Error fetching url %s", url)
        flash("Error fetching url")
        return redirect(url_for('start') + '?active=url')


def pastein():
    text = request.form.get('pastein')
    # email = request.args.get['email']
    if not text:
        flash('Text Required in Pastefield')
        return redirect(url_for('start') + '?active=pastein')
    (fd, savepath, relpath) = mksavefile("txt")
    logger.info("Saving paste to %r", relpath)
    with os.fdopen(fd, 'wb') as fh:
        fh.write(text)
    flash('redirecting ' + relpath)
    return redirect(url_for('run_frame', path=relpath))


def search():
    search = request.form.get('entrez_search_term')
    if not search:
        flash('Accession Number Required')
        return redirect(url_for('start') + '?active=search')
    self.active = 'entrez_search_term'
    self.session = entrez.CachedEntrezSession(library_root())

    self.session.search(cleaned_data['entrez_search_term'])
    logger.debug(
        "Search finished, found %d matches",
        self.session.result_count)
    if self.session.result_count == 1:
        cleaned_data['path'] = getrelpath(self.session.fetch())
    elif self.session.result_count > 1:
        raise Form.ValidationError(
            "Too many results (%d) found, need 1."
            " Try refining the search or searching for a RefSeq id."
            % (self.session.result_count))
    else:
        raise Form.ValidationError("No results found.")


def view():
    action = request.form.get('action')

    if request.method == 'POST':

        logger.info("Handling post action %r", action)
        if action == 'pastein':
            return pastein()
        elif action == 'fetchurl':
            return fetchurl()
        elif action == 'upload':
            return file_upload()
    #     if form.is_valid():

    #         logger.info("Form is valid; action is %r", action)
    #         kwargs = startform.cleaned_data
    #         path = kwargs.pop('path')
    #         if action == 'run':
    #             return redirect(url_for(action, path=path) + kwargs)
    #         else:
    #             logger.error("Unknown action %r", action)
    #             flash("Unknown action.")
    # else:
    #     email = startform['email'].data or request.args.get('email')
    return flask.render_template(
        'start.html', **{
            'action': action
        })


def re_search():
    if request.REQUEST.get('entrez_search_term'):
        return flask.render_template(
            'start.html')
    else:
        return view(request)


def efetch(id):
    logger.info("Asked to fetch Id: %s", id)
    searchsession = entrez.EntrezSession(library_root())
    abspath = searchsession.fetch_id(id)
    path = getrelpath(abspath)
    try:
        parsing.initial(abspath)
    except:
        flash("There was a problem loading file '%s', "
              "please try again or try a different record."
              % path)
        return re_search(request)

    path = getrelpath(abspath)
    remaining_args = dict(request.REQUEST.items())
    action = remaining_args.pop('action', 'run')

    if action in ['run', 'config']:
        return redirect(
            url_for(action, path=path) + remaining_args)
    else:
        logger.error("Unknown action %r", action)
        flash("Unknown action.")
        return re_search(request)
