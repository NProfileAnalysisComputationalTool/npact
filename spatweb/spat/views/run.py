# Create your views here.
import os.path, logging, tempfile

from django.conf import settings
from django.http import HttpResponse,HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.template import RequestContext
from django import forms
from django.contrib import messages


from __init__ import session_key, is_clean_path, getabspath, getrelpath
from spat.middleware import RedirectException
from spat import helpers
#from spat.helpers import add_help_text

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
    significance=forms.ChoiceField(choices=prepare.significance_levels)
    start_page=forms.IntegerField(required=False)
    end_page=forms.IntegerField(required=False)

def get_display_items(request, config):
    yield ('Filename', config['basename'])
    for key in ['date','length','description'] :
        if config.get(key) :
            yield key, config.get(key)


def run_it(request, path, form, config):
    logger.info("Got clean post, running.")
  
    config.update(form.cleaned_data)

    gbp = main.GenBankProcessor(getabspath(path), config=config)
    psname = gbp.run_Allplots()
    logger.debug("Got back ps file: %r", psname)
    
    psname = getrelpath(psname)
    logger.debug("relpath: %r",psname)
    raise RedirectException(reverse('results', args=[psname]))




def view(request, path):
    if not is_clean_path(path) :
        messages.error(request, "Path contained illegal characters, please upload a file or go to the library and select one.")
        return HttpResponseRedirect(reverse('spat.views.start.view'))

    form = None
    config = None
    try:
        config = prepare.default_config(getabspath(path))
    except:
        messages.error(request,"There was a problem loading file '%s', please try again or try a different record." % path)
        return HttpResponseRedirect(reverse('spat.views.start.view'))

    if request.method == 'POST' :
        form= RunForm(request.POST)
        if form.is_valid() :
            run_it(request, path, form, config)
    else :
        form = RunForm(initial=config)
    helpers.add_help_text(form,prepare.CONFIG_HELP_TEXT)

    return render_to_response('run.html',{'form':form, 'parse_data':config,
                                          'def_list_items': get_display_items(request,config)},
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
