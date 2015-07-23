import flask
import os.path
import tempfile
import urllib2
from npactflask import app
from npactflask.views import getrelpath, is_clean_path, library_root
from pynpact import util, entrez, parsing
from flask import url_for, redirect, request, flash


logger = app.logger


def mksavefile(prefix):
    """Wrapper around creating the file to save uploaded files in.

    Returns the (fd, abspath, relpath)
    """
    # we're using tempfile to ensure the file we get is unique and
    # aren't overwriting another.
    fd, abspath = tempfile.mkstemp(dir=app.config['UPLOADS'], prefix=prefix)
    relpath = getrelpath(abspath)
    return (fd, abspath, relpath)


def file_upload():
    fileup = request.files.get('file')
    args = dict(request.args.items())
    if not fileup:
        flash('Must upload a File')
        return redirect(url_for('start', **args))
    logger.info("Checking Uploaded file %s", fileup)
    if not is_clean_path(fileup.name):
        flash("Illegal filename")
        return redirect(url_for('start', **args))
    fd, savepath, relpath = mksavefile("up-%s-" % fileup.filename)
    logger.info("Saving uploaded file to %r", savepath)
    fileup.save(savepath)
    return redirect(url_for('run_frame', path=relpath))


def fetchurl():
    args = dict(request.values.items())
    url = args.pop('fetchurl', None)
    if not url:
        flash('URL Required in URL Field')
        return redirect(url_for('start', **args))
    logger.debug("Going to pull from %r", url)
    pull_req = urllib2.Request(url)
    if pull_req.get_type == "file":
        flash("Illegal URL")
        return redirect(url_for('start', **args))
    try:
        fh = urllib2.urlopen(pull_req)
        fd, savepath, relpath = mksavefile("url")
        util.stream_to_file(fh, savepath)
        return redirect(url_for('run_frame', path=relpath))
    except:
        logger.exception("Error fetching url %s", url)
        flash("Error fetching url")
        return redirect(url_for('start', **args))


def pastein():
    args = dict(request.values.items())
    text = args.pop('pastein', None)

    if not text:
        flash('Text Required in Pastefield')
        return redirect(url_for('start', **args))
    (fd, savepath, relpath) = mksavefile("txt")
    logger.info("Saving paste to %r", relpath)
    with os.fdopen(fd, 'wb') as fh:
        fh.write(text)
    flash('redirecting ' + relpath)
    return redirect(url_for('run_frame', path=relpath, **args))


def search():
    args = dict(request.values.items())
    search = args.pop('search', None)
    if not search:
        flash('Accession Number Required')
        return redirect(url_for('start', **args))
    sess = entrez.CachedEntrezSession(library_root())

    sess.search(search)
    logger.debug(
        "Search finished, found %d matches",
        sess.result_count)
    if sess.result_count == 1:
        relpath = getrelpath(sess.fetch())
        return redirect(url_for('run_frame', path=relpath, **args))
    elif sess.result_count > 1:
        flash(
            "Too many results (%d) found, need 1."
            " Try refining the search or searching for a RefSeq id."
            % (sess.result_count))
        return flask.render_template('start.html', sess=sess, search=search)
    else:
        flash("No results found.")
        return redirect(url_for('start', **args))


def view():
    active = request.form.get('active')
    if request.method == 'POST':

        logger.info("Handling post active %r", active)
        if active == 'pastein':
            return pastein()
        elif active == 'fetchurl':
            return fetchurl()
        elif active == 'upload':
            return file_upload()
        elif active == 'search':
            return search()
    else:
        # not a post
        args = dict(request.args.items())
        return flask.render_template('start.html', **args)


def re_search():
    if request.args.get('entrez_search_term'):
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
        return re_search()

    args = dict(request.args.items())

    return redirect(url_for('run_frame', path=path, **args))
