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


def get_ti(**kwargs) :
    return forms.TextInput(attrs=kwargs)

class RunForm(forms.Form) :
    first_page_title = forms.CharField(widget=get_ti(size=40))
    following_page_title = forms.CharField(required=False,widget=get_ti(size=40))
    length=forms.IntegerField(required=True, min_value=0,
                              widget=get_ti(size=8))

def prefill_form(request, path, parse_data) :

    title = parse_data.get('description') or parse_data.get('basename')
    form = RunForm(initial={'first_page_title': title,
                    'following_page_title': title,
                    'length': parse_data.get('length')
                    })

    return form

def get_display_items(request,parse_data) :
    yield ('Filename', parse_data['basename'])
    for key in ['date','length','description'] :
        if parse_data.get(key) :
            yield key, parse_data.get(key)


def run_it(request, path, form,parse_data) :
    config = prepare.default_config(parse_data)
    for k in ['first_page_title', 'following_page_title'] :
        if form.cleaned_data.get(k) :
            config[k] = form.cleaned_data[k]

    gbp = main.GenBankProcessor(os.path.join(settings.MEDIA_ROOT, path), config=config)
    psname = gbp.run_Allplots()
    logger.debug("Got back ps file: %r", psname)
    psname = os.path.relpath(psname, settings.MEDIA_ROOT)
    logger.debug("relpath: %r",psname)
    raise RedirectException(reverse('results', args=[psname]))




def view(request, path):
    if not is_clean_path(path) :
        messages.error(request, "Path contained illegal characters, please upload a file or go to the library and select one.")
        return HttpResponseRedirect(reverse('spat.views.start.view'))

    form = None
    abs_path = os.path.join(settings.MEDIA_ROOT, path)
    logger.debug("try_parseing %r", abs_path)
    parse_data = prepare.try_parse(abs_path)
    if not parse_data:
        messages.error(request,"There was a problem loading file '%s', please try again or try a different record." % path)
        return HttpResponseRedirect(reverse('spat.views.start.view'))

    if request.method == 'POST' :
        form= RunForm(request.POST)
        if form.is_valid() :
            run_it(request, path, form, parse_data)
    else :
        form = prefill_form(request, path, parse_data)

    return render_to_response('run.html',{'form':form, 'parse_data':parse_data,
                                          'def_list_items': get_display_items(request,parse_data)},
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
