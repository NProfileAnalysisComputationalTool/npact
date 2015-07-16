import flask
import logging
import os.path
import tempfile
import urllib2
from flask_wtf import Form, file, html5
from wtforms import TextField, TextAreaField
from npactflask.views import settings, getrelpath
from npactflask.views import is_clean_path, library_root
from pynpact import util, entrez, parsing
from flask import url_for, redirect, request, flash


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
    return TextField(attrs={'size': size})


class StartForm(Form):
    active = None
    accordion_fields = ['file_upload', 'url', 'pastein', 'entrez_search_term']
    file_upload = file.FileField()
    url = html5.URLField(label="From URL",
                         widget=get_ti(50))
    pastein = TextField(
        label="Paste in as text",
        widget=TextAreaField(attrs={'rows': 3, 'cols': ""}))
    entrez_search_term = TextField(
        label="Accession Number")
    email = html5.EmailField()

    def __init__(self, *args, **kwargs):
        super(StartForm, self).__init__(*args, **kwargs)

    def clean_file_upload(self):
        cleaned_data = self.cleaned_data

        if cleaned_data.get('file_upload'):
            self.active = 'file_upload'
            fu = cleaned_data.get('file_upload')
            logger.info("Checking Uploaded file %s", fu)
            if not is_clean_path(fu.name):
                raise Form.ValidationError("Illegal filename")

            fd, savepath, relpath = mksavefile("up-%s-" % fu.name)
            logger.info("Saving uploaded file to %r", relpath)
            with os.fdopen(fd, 'wb') as fh:
                for chunk in fu.chunks():
                    fh.write(chunk)

            cleaned_data['path'] = relpath

    def clean_url(self):
        cleaned_data = self.cleaned_data
        if cleaned_data.get('url'):
            self.active = 'url'
            url = cleaned_data.get('url')
            logger.debug("Going to pull from %r", url)
            pull_req = urllib2.Request(url)
            if pull_req.get_type == "file":
                raise Form.ValidationError("Illegal URL")
            try:
                fh = urllib2.urlopen(pull_req)
                fd, savepath, relpath = mksavefile("url")
                util.stream_to_file(fh, savepath)
                cleaned_data['path'] = relpath
            except:
                logger.exception("Error fetching url %s", url)
                raise Form.ValidationError("Error fetching url")

    def clean_pastein(self):
        cleaned_data = self.cleaned_data
        if cleaned_data.get('pastein'):
            self.active = 'pastein'
            text = cleaned_data.get('pastein')
            (fd, savepath, relpath) = mksavefile("txt")
            logger.info("Saving paste to %r", relpath)
            cleaned_data['path'] = relpath
            with os.fdopen(fd, 'wb') as fh:
                fh.write(text)

    def clean_entrez_search_term(self):
        cleaned_data = self.cleaned_data
        if cleaned_data.get('entrez_search_term'):
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

    def clean(self):
        cleaned_data = self.cleaned_data

        if not cleaned_data.get('path'):
            raise Form.ValidationError("No valid GenBank file submitted.")
        return cleaned_data


def view():
    action = None
    if request.method == 'POST':
        path = None
        startform = StartForm(request.POST, request.FILES)

        action = request.POST.get('action')
        if startform.is_valid():

            logger.info("Form is valid; action is %r", action)
            kwargs = startform.cleaned_data
            path = kwargs.pop('path')
            if action == 'run':
                return redirect(url_for(action, path=path) + kwargs)
            else:
                logger.error("Unknown action %r", action)
                flash("Unknown action.")
    else:
        startform = StartForm()
        email = startform['email'].data or request.args.get('email')
    return flask.render_template(
        'start.html', **{
            'form': startform,
            'action': action,
            'email': email
        })


def re_search():
    if request.REQUEST.get('entrez_search_term'):
        return flask.render_template(
            'start.html', **{'form': StartForm(request.REQUEST)})
    else:
        return view(request)


def efetch(id):
    logger.info("Asked to fetch Id: %s", id)
    session = entrez.EntrezSession(library_root())
    abspath = session.fetch_id(id)
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
