# Create your views here.
import logging
import os.path
import tempfile
import urllib2

from pynpact import util, entrez, prepare

from django import forms
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect
from django.shortcuts import render_to_response
from django.template import RequestContext


from spatweb import is_clean_path, getabspath, getrelpath, library_root
from spatweb.middleware import RedirectException


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



class StartForm(forms.Form):
    active = None
    file_upload = forms.FileField(required=False)
    url = forms.URLField(label="From URL",required=False)
    pastein = forms.CharField(label="Paste in as text",
                              widget=forms.Textarea(attrs={'rows':3, 'cols':""}),
                              required=False)
    entrez_search_term = forms.CharField(label="Search Entrez for Term", required=False)

    def __init__(self,*args,**kwargs):
        super(StartForm,self).__init__(*args,**kwargs)


    def clean_file_upload(self):
        cleaned_data = self.cleaned_data

        if cleaned_data.get('file_upload'):
            self.active = 'file_upload'
            fu = cleaned_data.get('file_upload')
            logger.info("Checking Uploaded file %s", fu)
            if not is_clean_path(fu.name):
                raise forms.ValidationError("Illegal filename")

            fd,savepath,relpath = mksavefile( "up-%s-" % fu.name)
            logger.info("Saving uploaded file to %r", relpath)
            with os.fdopen(fd, 'wb') as fh :
                for chunk in fu.chunks() :
                    fh.write(chunk)

            cleaned_data['path']=relpath

    def clean_url(self):
        cleaned_data = self.cleaned_data
        if cleaned_data.get('url'):
            self.active='url'
            url = cleaned_data.get('url')
            logger.debug("Going to pull from %r", url)
            pull_req = urllib2.Request(url)
            if pull_req.get_type == "file":
                raise forms.ValidationError("Illegal URL")
            try:
                fh = urllib2.urlopen(pull_req)
                fd,savepath,relpath = mksavefile("url")
                util.stream_to_file(fh,savepath)
                cleaned_data['path'] = relpath
            except:
                logger.exception("Error fetching url %s", url)
                raise forms.ValidationError("Error fetching url")


    def clean_pastein(self):
        cleaned_data = self.cleaned_data
        if cleaned_data.get('pastein'):
            self.active = 'pastein'
            text = cleaned_data.get('pastein')
            (fd,savepath,relpath) = mksavefile("txt")
            logger.info("Saving paste to %r", relpath)
            cleaned_data['path'] = relpath
            with os.fdopen(fd,'wb') as fh :
                fh.write(text)

    def clean_entrez_search_term(self):
        cleaned_data = self.cleaned_data
        if cleaned_data.get('entrez_search_term'):
            self.active='entrez_search_term'
            self.session = entrez.CachedEntrezSession(library_root())

            self.session.search(cleaned_data['entrez_search_term'])
            logger.debug("Search finished, found %d matches", self.session.result_count)
            if self.session.result_count == 1:
                cleaned_data['path'] = getrelpath(self.session.fetch())
            elif self.session.result_count > 1:
                raise forms.ValidationError("Too many results (%d) found, need 1."
                                            " Try refining the search or searching for a RefSeq id."
                                            % (self.session.result_count))
            else:
                raise forms.ValidationError("No results found.")

    def clean(self):
        cleaned_data = self.cleaned_data

        if not cleaned_data.get('path'):
            raise forms.ValidationError("No valid GenBank file submitted.")
        return cleaned_data


def view(req) :
    form = None
    if req.method == 'POST' :
        path=None
        startform = StartForm(req.POST, req.FILES)

        if startform.is_valid():
            logger.info("Form is valid; action is %r", req.POST.get('action'))
            path = startform.cleaned_data['path']
            return HttpResponseRedirect(reverse('run',args=[path]))
    else:
        startform = StartForm()

    return render_to_response('start.html', {'form': startform},
                              context_instance=RequestContext(req))


def re_search(req):
    if req.REQUEST.get('entrez_search_term'):
        return render_to_response('start.html', {'form': StartForm(req.REQUEST)},
                                  context_instance=RequestContext(req))
    else:
        return view(req)



def efetch(req, id):
    logger.info("Asked to fetch Id: %s", id)
    session = entrez.EntrezSession(library_root())
    abspath = session.fetch_id(id)

    try:
        prepare.try_parse(abspath)
    except prepare.InvalidGBKException, e:
        messages.error(req, str(e))
        return re_search(req)
    except:
        messages.error(request,
                       "There was a problem loading file '%s', "
                       "please try again or try a different record."
                       % path)
        return re_search(req)

    return HttpResponseRedirect(reverse('run',args=[getrelpath(abspath)]))
