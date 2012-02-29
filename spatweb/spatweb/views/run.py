# Create your views here.
import logging
import os.path
import json

from django import forms
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from django.utils.http import urlencode
from django.template import RequestContext
from pynpact import prepare, main, util
from pynpact.softtimeout import Timeout

from spatweb import assert_clean_path, getabspath, getrelpath
from spatweb import helpers
from spatweb.middleware import RedirectException

#from spatweb.helpers import add_help_text

# Get an instance of a logger
logger = logging.getLogger(__name__)



def get_reconfigure_url(request, path=None):
    if path is None:
        path = request.GET.get('path')
    if path:
        return reverse('config', args=[path]) + "?" + urlencode(request.GET)
    else:
        return None

def get_raw_url(request, path):
    #return request.build_absolute_uri(reverse('raw', path))
    return reverse('raw', args=[path])


def get_ti(size):
    return forms.TextInput(attrs={'size':size})

class ConfigForm(forms.Form):
    first_page_title = forms.CharField(widget=get_ti(40))
    following_page_title = forms.CharField(required=False, widget=get_ti(40))
    length=forms.IntegerField(required=True, min_value=0,
                              widget=get_ti(8))
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
    for key in ['date','length','description']:
        if config.get(key):
            yield key, config.get(key)

        
def config(request, path):
    assert_clean_path(path, request)
    config = build_config(path, request)

    form = None
    if request.method == 'POST':
        form = ConfigForm(request.POST)
        if form.is_valid():
            logger.info("Got clean post, running.")
            url = reverse('run', args=[path]) + "?" + urlencode(form.cleaned_data)
            return HttpResponseRedirect(url)
    else:
        form = ConfigForm(initial=config)

    helpers.add_help_text(form, prepare.CONFIG_HELP_TEXT)

    return render_to_response('config.html',
                              {'form':form, 'parse_data':config,
                               'def_list_items': get_display_items(request,config)},
                               context_instance=RequestContext(request))

def build_config(path, request):
    assert_clean_path(path, request)

    try:
        config = prepare.default_config(getabspath(path))
    except prepare.InvalidGBKException, e:
        messages.error(request, str(e))
        raise RedirectException(reverse('start'))
    except:
        logger.exception("Error parsing gbk: %r", getabspath(path))
        messages.error(request,
                       "There was a problem loading file '%s', "
                       "please try again or try a different record." % path)
        raise RedirectException(reverse('start'))

    cf = ConfigForm(request.REQUEST)
    for f in cf.visible_fields():
        try:
            v=f.field.clean(f.field.to_python(f.data))
            if v:
                logger.debug("Including %r:%r from request.", f.name, v)
                config[f.name] = v
        except forms.ValidationError, ve:
            pass
        except:
            logger.exception("Error with %r", f.name)
    if cf.is_valid():
        logger.debug('updating with %r', cf.cleaned_data)
        config.update(cf.cleaned_data)

    return config
    
def encode_config(config, **urlconf):
    for k in ConfigForm().fields.keys():
        v = config.get(k, None)
        if v:
            urlconf[k] = v
    return "?" + urlencode(urlconf)


def run_frame(request, path):
    config = build_config(path, request)
    request.session[path] = config
    full_path = reverse('process',args=[path])
    return render_to_response('processing.html',
                              {'path': full_path},
                              context_instance=RequestContext(request))


def run_step(request, path):
    """Invoked via ajax, runs part of the process with a softtimeout until finished."""
    assert_clean_path(path, request)
    try: 
        config = request.session.get(path)
        #the frame is supposed to ensure this is in session.
        if not config:
            return HttpResponse('Session Timeout, please try again.', status=500)

        gbp = main.GenBankProcessor(getabspath(path), config=config, timeout=4)
        nextstep = None
        try:
            pspath = gbp.process()
            logger.debug("Finished processing.")
            pspath = getrelpath(pspath)
            #url = reverse('results', args=[psname]) + encode_config(config, path=path)
            nextstep = {'next':'results', 
                        'download_url': get_raw_url(request, pspath),
                        'reconfigure_url': reverse('config', args=[path]) + encode_config(config)}
        except Timeout, pt:
            nextstep = {'next':'process', 'pt': vars(pt)}
        return HttpResponse(json.dumps(nextstep))
    except:
        logger.exception("Error in run_step")
        return HttpResponse('ERROR', status=500)


def results(request, path):
    """Serve a results page."""
    assert_clean_path(path, request)

    download_link = None
    try:
        getabspath(path)
        download_link=get_raw_url(request, path)
    except IOError:
        messages.error(request,
                       "We're sorry but that file no longer exists. We "
                       "delete old results periodically to save space on"
                       " the server. Please try running the analysis again.")

    return_url = get_return_url(request)

    return render_to_response('results.html',
                              {'download_link': download_link,
                               'return_url': return_url},
                              context_instance=RequestContext(request))
