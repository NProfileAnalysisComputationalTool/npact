# Create your views here.
import os.path, logging, tempfile

from django.conf import settings
from django.http import HttpResponse,HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.template import RequestContext
from django import forms
from django.contrib import messages


from __init__ import session_key, is_clean_path
from spat.middleware import RedirectException

from pynpact import prepare, main


# Get an instance of a logger
logger = logging.getLogger(__name__)


class RunForm(forms.Form) :
    first_page_title = forms.CharField()
    following_page_title = forms.CharField(required=False)


def prefill_form(request, path, data) :

    title = data.get('description') or data.get('basename')
    form = RunForm({'first_page_title': title,
                    'following_page_title': title})

    return form

def get_display_items(request,data) :
    yield ('Filename', data['basename'])
    for key in ['date','length','description'] :
        if data.get(key) :
            yield key, data.get(key)


def run_it(request, path, form) :
    config = {}
    for k in ['first_page_title', 'following_page_title'] :
        if form.cleaned_data.get(k) :
            config[k] = form.cleaned_data[k]

    gbp = main.GenBankProcessor(os.path.join(settings.MEDIA_ROOT, path), config=config)
    psname = gbp.run_Allplots()
    psname = os.path.relpath(psname, settings.MEDIA_ROOT)
    raise RedirectException(reverse('results', args=[psname]))




def view(request, path):
    if not is_clean_path(path) :
        messages.error(request, "Path contained illegal characters, please upload a file or go to the library and select one.")
        return HttpResponseRedirect(reverse('spat.views.start.view'))

    form = None
    data = prepare.try_parse(os.path.join(settings.MEDIA_ROOT, path))
    if not data:
        messages.error(request,"There was a problem loading file '%s', please try again or try a different record." % path)
        return HttpResponseRedirect(reverse('spat.views.start.view'))

    if request.method == 'POST' :
        form= RunForm(request.POST)
        if form.is_valid() :
            run_it(request, path, form)
    else :
        form = prefill_form(request, path, data)

    return render_to_response('run.html',{'form':form, 'data':data,
                                          'def_list_items': get_display_items(request,data)},
                               context_instance=RequestContext(request))


def view_none(request) :
    messages.error(request, "No genome source selected, please upload one, or go to the library and select one.")
    return HttpResponseRedirect(reverse('spat.views.start.view'))





def results(request,path) :
    """Serve either a results page, or if the 'raw' querystring parameter is set to anything then return the ps file directly."""

    if not is_clean_path(path) :
        messages.error(request, "Path contained illegal characters, please upload a file or go to the library and select one.")
        return HttpResponseRedirect(reverse('start'))

    download_link=request.build_absolute_uri(reverse('raw', kwargs={'path':path}))
    return render_to_response('results.html',
                              {'download_link': download_link},
                              context_instance=RequestContext(request))
