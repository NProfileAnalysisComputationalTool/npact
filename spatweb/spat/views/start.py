# Create your views here.
import os.path, logging, tempfile

from django.conf import settings
from django.http import HttpResponse,HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.template import RequestContext
from django import forms
from django.contrib import messages

from __init__ import session_key

# Get an instance of a logger
logger = logging.getLogger(__name__)


class StartForm(forms.Form):
    file_upload = forms.FileField(required=False)
    url = forms.URLField(label="From URL", required=False)
    pastein = forms.CharField(label="Paste it in",widget=forms.Textarea,required=False)


def saveFile(request, fu) :
    savepath = os.path.join(settings.MEDIA_ROOT,fu.name)
    logger.info("Saving uploaded file to %r", savepath)
    #TODO: this should probably be sanitized, and make sure we aren't stepping on someone else.
    with open(savepath, 'wb') as fh :
        for chunk in fu.chunks() :
            fh.write(chunk)

    return savepath

def pullFromUrl(request, url) :
    raise NotImplemented("Pulling a file from a url is not yet implemented.")


def saveToFile(request, text) :
    (fd,savepath) = tempfile.mkstemp(dir=settings.MEDIA_ROOT, prefix="txt", suffix=".gbk")
    logger.info("Saving paste to %r, %r", savepath, __name__)
    with os.fdopen(fd,'wb') as fh :
        fh.write(text)
    return savepath




def view(request) :
    form = None
    path = None
    if request.method == 'POST' :
        form=StartForm(request.POST, request.FILES)
        if form.is_valid() :
            if request.FILES.get('file_upload') :
                path = saveFile(request, request.FILES.get('file_upload'))
            elif form.cleaned_data.get('url') :
                path = pullFromUrl(request, form.cleaned_data.get('url'))
            elif form.cleaned_data.get('pastein') :
                path = saveToFile(request, form.cleaned_data.get('pastein'))

            if path :
                request.session[session_key("input_path")] = path
                return HttpResponseRedirect(reverse('spat.views.run.view'))
            else :
                messages.error(request, "You need to supply a genome source.")



    return render_to_response('start.html',{'form':form or StartForm()},
                                  context_instance=RequestContext(request))
