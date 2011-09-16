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


def mksavefile(prefix):
    """Wrapper around creating the file to save uploaded files in.

  Returns the (fd,path) tuple (similar to tempfile.mkstemp)
"""
    #we're using tempfile to ensure the file we get is unique and
    #aren't overwriting another.
    fd,abspath = tempfile.mkstemp(dir=settings.MEDIA_ROOT, prefix=prefix, suffix=".gbk")
    relpath = getrelpath(abspath)
    return (fd,abspath,relpath)



            # elif form.cleaned_data.get('url') :
            #     path = pullFromUrl(req, form.cleaned_data.get('url'))
            # elif form.cleaned_data.get('pastein') :
            #     path = saveToFile(req, form.cleaned_data.get('pastein'))


class UploadForm(forms.Form):
    type="UploadForm"
    title="Upload a File"
    file_upload = forms.FileField(required=True)

    def clean(self):
        cleaned_data = self.cleaned_data
        fu = cleaned_data['file_upload']
        logger.info("Checking Uploaded file %s", fu)
        if not is_clean_path(fu.name):
            raise forms.ValidationError("Illegal filename")

        fd,savepath,relpath = mksavefile( "up-%s-" % fu.name)
        logger.info("Saving uploaded file to %r", relpath)
        with os.fdopen(fd, 'wb') as fh :
            for chunk in fu.chunks() :
                fh.write(chunk)

        cleaned_data['path']=relpath
        return cleaned_data

class UrlForm(forms.Form):
    type="UrlForm"
    title="From URL"
    url = forms.URLField(label="From URL", required=True)

    def clean(self):
        cleaned_data = self.cleaned_data
        url = cleaned_data['url']
        logger.debug("Going to pull from %r", url)
        pull_req = urllib2.Request(url)
        if pull_req.get_type == "file":
            raise forms.ValidationError("Illegal URL")
        try:
            fh = urllib2.urlopen(pull_req)
            fd,savepath,relpath = mksavefile("url")
            pynpact.util.stream_to_file(fh,savepath)
            cleaned_data['path'] = relpath
        except:
            logger.exception("Error fetching url %s", url)
            raise forms.ValidationError("Error fetching url")

        return cleaned_data

class PasteForm(forms.Form):
    type="PasteForm"
    title="Paste as Text"
    pastein = forms.CharField(label="Paste in as text",
                              widget=forms.Textarea(attrs={'rows':3, 'cols':""}),
                              required=True)
    def clean(self):
        cleaned_data = self.cleaned_data
        text = cleaned_data['pastein']
        (fd,savepath,relpath) = mksavefile("txt")
        logger.info("Saving paste to %r", relpath)
        cleaned_data['path'] = relpath
        with os.fdopen(fd,'wb') as fh :
            fh.write(text)
        return cleaned_data

class EntrezSearchForm(forms.Form):
    type="EntrezSearchForm"
    title="Search Entrez"
    term = forms.CharField(label="Search Term")
    refine = forms.BooleanField("Refine previous search",
                                initial=True,
                                help_text="Leave checked to search inside the previous results.")
    def clean(self):
        pass

def build_forms(req):
    types = [UploadForm,UrlForm,PasteForm,EntrezSearchForm]
    if req.method == 'POST':
        forms = []
        for type in types:
            if type.__name__ == req.POST['type']:
                forms.insert(0, type(req.POST,req.FILES))
            else:
                forms.append(type())
        return forms
    else:
        return [type() for type in types]


def view(req) :
    form = None
    if req.method == 'POST' :
        path=None
        forms=build_forms(req)
        submitted = forms[0]
        if submitted.is_valid():
            logger.info("Form is valid")
            path = submitted.cleaned_data['path']
            raise RedirectException(reverse('run',args=[path]))
    else:
        forms = build_forms(req)

    return render_to_response('start.html',{'forms':forms},
                                  context_instance=RequestContext(req))
