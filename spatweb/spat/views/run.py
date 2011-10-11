# Create your views here.
import logging
import os.path

from django import forms
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect
from django.shortcuts import render_to_response
from django.utils.http import urlencode
from django.template import RequestContext
from pynpact import prepare, main, util
from spat import is_clean_path, getabspath, getrelpath
from spat import helpers
from spat.middleware import RedirectException
from spat.views import get_return_url

#from spat.helpers import add_help_text


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

    def clean(self):
        cleaned_data = self.cleaned_data
        start = cleaned_data.get('start_page')
        end = cleaned_data.get('end_page')
        if start and end and start > end:
            raise forms.ValidationError("End page must be greater than or equal to start page.")
        return cleaned_data

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
    url = reverse('results', args=[psname])
    urlconf = util.reducedict(config, form.fields.keys())
    urlconf['path'] = path
    
    url += "?" + urlencode(urlconf)
    raise RedirectException(url)




def view(request, path):
    if not is_clean_path(path) :
        messages.error(request, "Path contained illegal characters, please upload a file or go to the library and select one.")
        return HttpResponseRedirect(reverse('spat.views.start.view'))

    form = None
    config = None
    try:
        config = prepare.default_config(getabspath(path))
    except prepare.InvalidGBKException, e:
        messages.error(request,str(e))
        return HttpResponseRedirect(reverse('spat.views.start.view'))
    except:
        messages.error(request,"There was a problem loading file '%s', please try again or try a different record." % path)
        return HttpResponseRedirect(reverse('spat.views.start.view'))

    if request.method == 'POST' :
        form = RunForm(request.POST)
        if form.is_valid():
            run_it(request, path, form, config)
    else:
        #can't use dict.update here as the multi-value dict gives back
        #an array in that case
        for k in request.GET:
            logger.info('k:%s, v:%s', k, request.GET[k])
            config[k]= request.GET[k]
        form  = RunForm(initial=config)

    helpers.add_help_text(form,prepare.CONFIG_HELP_TEXT)

    return render_to_response('run.html',{'form':form, 'parse_data':config,
                                          'def_list_items': get_display_items(request,config)},
                               context_instance=RequestContext(request))


def view_none(request) :
    messages.error(request, "No genome source selected, please upload one, or go to the library and select one.")
    return HttpResponseRedirect(reverse('spat.views.start.view'))





def results(request, path):
    """Serve either a results page, or if the 'raw' querystring
    parameter is set to anything then return the ps file directly."""

    if not is_clean_path(path) :
        messages.error(request, "Path contained illegal characters, please upload a file or go to the library and select one.")
        return HttpResponseRedirect(reverse('start'))

    download_link = None
    try:
        getabspath(path)
        download_link=request.build_absolute_uri(reverse('raw', kwargs={'path':path}))
    except IOError:
        messages.error(request, "We're sorry but that file no longer exists. We expire old results periodically to save space on the server. Please try running the analysis again.")

    return_url = get_return_url(request)

    return render_to_response('results.html',
                              {'download_link': download_link,
                               'return_url': return_url},
                              context_instance=RequestContext(request))
