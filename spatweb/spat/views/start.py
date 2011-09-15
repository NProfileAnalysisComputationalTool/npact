# Create your views here.
import logging
import os.path
import tempfile
import urllib2

import pynpact.util

from django import forms
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect
from django.shortcuts import render_to_response
from django.template import RequestContext

from spat.middleware import RedirectException
from spat.views import is_clean_path, getabspath, getrelpath


# Get an instance of a logger
logger = logging.getLogger(__name__)

class UploadForm(forms.Form):
    file_upload = forms.FileField(required=False)

class UrlForm(forms.Form):
    url = forms.URLField(label="From URL", required=False)

class PasteForm(forms.Form):
    pastein = forms.CharField(label="Paste it in",widget=forms.Textarea,required=False)

class EntrezSearchForm(forms.Form):
    term = forms.CharField(label="Search Term")
    refine = forms.BooleanField("Refine previous search",
                                initial=True,
                                help_text="Leave checked to search inside the previous results.")

def mksavefile(req, prefix):
    """Wrapper around creating the file to save uploaded files in.

  Returns the (fd,path) tuple (similar to tempfile.mkstemp)
"""
    #we're using tempfile to ensure the file we get is unique and
    #aren't overwriting another.
    fd,abspath = tempfile.mkstemp(dir=settings.MEDIA_ROOT, prefix=prefix, suffix=".gbk")
    relpath = getrelpath(abspath)
    return (fd,abspath,relpath)

def saveUploadFile(req, fu) :
    if not is_clean_path(fu.name) :
        messages.error(req,"Illegal filename")
        raise IOError("Illegal filename")

    fd,savepath,relpath = mksavefile(req, "up-" + fu.name)
    logger.info("Saving uploaded file to %r", relpath)
    with os.fdopen(fd, 'wb') as fh :
        for chunk in fu.chunks() :
            fh.write(chunk)

    return relpath

def pullFromUrl(req, url) :
    logger.debug("Going to pull from %r", url)
    pull_req = urllib2.Request(url)
    if pull_req.get_type == "file":
        messages.error(req,"Illegal URL.")
        raise IOError("Illegal URL")
    try:
        with urllib2.urlopen(pull_req) as fh:
            fd,savepath,relpath = mksavefile(req,"url")
            pynpact.util.stream_to_file(fh,savepath)
    except:
        messages.error(req,"Error opening url.")
        raise

def saveToFile(req, text) :
    (fd,savepath,relpath) = mksavefile(req,"txt")
    logger.info("Saving paste to %r", relpath)
    with os.fdopen(fd,'wb') as fh :
        fh.write(text)
    return relpath


def view_POST(req):
    path=None
    form=StartForm(req.POST, req.FILES)
    entrez_search = EntrezSearchForm(req.POST)
    if form.is_valid() :
        try:
            if req.FILES.get('file_upload') :
                path = saveUploadFile(req, req.FILES.get('file_upload'))
            elif form.cleaned_data.get('url') :
                path = pullFromUrl(req, form.cleaned_data.get('url'))
            elif form.cleaned_data.get('pastein') :
            path = saveToFile(req, form.cleaned_data.get('pastein'))

            if path :
                return RedirectException(reverse('run',args=[path]))
            else :
                messages.error(req, "You need to supply a genome source.")
        except Exception,e:
            logger.exception(e)
    elif entrez_search.is_valid() :
        
    


def view(req) :
    form = None
    if req.method == 'POST' :
        form = view_POST(req)

    return render_to_response('start.html',{'form':form or StartForm()},
                                  context_instance=RequestContext(req))
