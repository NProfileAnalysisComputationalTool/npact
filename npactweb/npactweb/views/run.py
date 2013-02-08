# Create your views here.
import logging
import os
import os.path
import json
from datetime import datetime

from django import forms
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from django.utils.http import urlencode
from django.template import RequestContext, Context
from django.template.loader import get_template

from pynpact import prepare, util
from pynpact import main
from pynpact.softtimeout import Timeout

from npactweb import assert_clean_path, getabspath, getrelpath
from npactweb.middleware import RedirectException

from taskqueue import client, NoSuchTaskError



# Get an instance of a logger
logger = logging.getLogger(__name__)




def get_raw_url(request, path):
    #return request.build_absolute_uri(reverse('raw', path))
    if path.startswith('/'):
        path = getrelpath(path)
    return reverse('raw', args=[path])


def get_ti(size):
    return forms.TextInput(attrs={'size':size})

class ConfigForm(forms.Form):
    first_page_title = forms.CharField(widget=get_ti(40))
    following_page_title = forms.CharField(required=False, widget=get_ti(40))
    length=forms.IntegerField(required=True, min_value=0,
                              widget=get_ti(8))
    nucleotides=forms.MultipleChoiceField(choices=[(i,i) for i in ['a','c','g','t']],
                                          widget=forms.CheckboxSelectMultiple())
    skip_prediction=forms.BooleanField(required=False)
    significance=forms.ChoiceField(choices=prepare.significance_levels, required=False,
                                   label="Prediction Significance")
    start_base=forms.IntegerField()
    end_base=forms.IntegerField()
    alternate_colors=forms.BooleanField(required=False, label="Alternative colors")
    email = forms.EmailField(required=False)

    def clean(self):
        cleaned_data = self.cleaned_data
        start = cleaned_data.get('start_base')
        end = cleaned_data.get('end_base')
        if start and end and start > end:
            raise forms.ValidationError("End page must be greater than or equal to start page.")
        return cleaned_data

def get_display_items(request, config):
    yield ('Filename', config['basename'])
    for key in ['date','length','description']:
        if config.get(key):
            yield key, config.get(key)


def config(request, path):
    "The config view. Renders the configform for the given gbkpath."
    assert_clean_path(path, request)
    config = build_config(path, request)

    form = None
    if request.method == 'POST':
        form = ConfigForm(request.POST)
        if form.is_valid():
            logger.info("Got clean post, running.")
            url = reverse('run', args=[path]) + "?" + urlencode(form.cleaned_data, True)
            return HttpResponseRedirect(url)
    else:
        form = ConfigForm(initial=config)

    for key,field in form.fields.items():
        if key in prepare.CONFIG_HELP_TEXT:
            field.help_text = prepare.CONFIG_HELP_TEXT[key]
        elif settings.DEBUG:
            logger.error("Help text missing for config form field: %r", key)

    return render_to_response('config.html',
                              {'form':form, 'parse_data':config,
                               'def_list_items': get_display_items(request,config)},
                               context_instance=RequestContext(request))

#Variables, mostly for debugging, that aren't exposed anywhere but
#that if present in the request should be added to the config
MAGIC_PARAMS = ['raiseerror', 'force']

def build_config(path, request):
    "Tries to build the config dictionary for the given path"
    assert_clean_path(path, request)

    try:
        config = prepare.default_config(getabspath(path))
    except prepare.InvalidGBKException, e:
        messages.error(request, str(e))
        raise RedirectException(reverse('start'))
    except:
        logger.exception("Error parsing gbk: %r", getabspath(path, raise_on_missing=False))
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

    for key in MAGIC_PARAMS:
        if key in request.GET:
            config[key] = request.GET[key]

    return config


def urlencode_config(config, exclude=None):
    """Encode a config dictionary suitably for passing in a url.

    This doesn't need to have anything in it but what can be entered
    in the ConfigForm; everything else can be recalculated based on
    that.
    """
    keys = ConfigForm().fields.keys() + MAGIC_PARAMS
    if exclude:
        keys = set(keys) - set(exclude)
    return urlencode(util.reducedict(config, keys), True)

def run_frame(request, path):
    """This is the main processing page. However, that page is just
    frame in which ajax requests spur the actual processing"""
    assert_clean_path(path, request)
    config = build_config(path, request)
    jobid = kickstart(request, path, config)
    email = request.GET.get('email')

    if email:
        results_link = reverse('run', args=[path]) + '?' + urlencode_config(config, exclude=['email'])
        results_link = request.build_absolute_uri(results_link)
        
        
        # the config dictionary at the end of the process is what is
        # returned by the above job, will be passed as the first
        # argument to the function called by client.after: send_email
        email_jobid = client.after(jobid, send_email, [results_link, email])
        logger.debug("Schedule email to %r with jobid: %s", email, email_jobid)


    
    #TODO: we could wait a tiny amount of time here to see if the job is already done.

    full_path = reverse('runstatus',args=[jobid])
    reconfigure_url = reverse('config', args=[path]) + "?" + urlencode_config(config)
    return render_to_response('processing.html',
                              {'statusurl': full_path,
                               'reconfigure_url': reconfigure_url,
                               'email': email,
                               },
                              context_instance=RequestContext(request))

def run_status(request, jobid):
    config = None
    result = {}
    status = 200

    try:
        ready = client.ready(jobid)
    except NoSuchTaskError:
        result['message'] = 'Unknown task identifier. Please retry.'
        status = 500
    except Exception,e:
        result['message'] = 'Fatal error.'
        status = 500
        logger.exception("Error getting job status. %r", jobid)

    if status != 200:
        #something has already gone wrong.
        return HttpResponse(json.dumps(result), status=status)

    if not ready:
        result['steps'] = client.log_tail(jobid, 1)
    else:
        try:
            config = client.result(jobid)
        except NoSuchTaskError:
            result['message'] = 'Unknown task identifier. Please retry.'
            status = 500
        except Exception,e:
            result['message'] = 'Fatal error.'
            status = 500
            logger.exception("Error getting job result. %r", jobid)

    if config:
        pdf_url = get_raw_url(request,
                              getrelpath(config.get('pdf_filename') or
                                         config.get('combined_ps_name')))
        result['next'] = 'results'
        result['download_url'] = pdf_url

        ago = config.get('acgt_gamma_output')
        if ago:
            result['files'] = [get_raw_url(request, os.path.join(ago, v))
                               for v in os.listdir(ago)]

    return HttpResponse(json.dumps(result), status=status)


### This isn't used anymore (in favor of runstatus). Not deleted yet
### in case we want to try to keep this branch as a fallback path in
### case the processing daemon isn't running.
# def run_step(request, path):
#     """Invoked via ajax, runs part of the process with a softtimeout until finished."""
#     assert_clean_path(path, request)
#     gbp,config,result = None,None,None
#     status=200

#     try:
#         config = request.session.get(path)
#         #the frame is supposed to ensure this is in session.
#         if not config:
#             return HttpResponse('Session Timeout, please try again.', status=500)
#         gbp = main.GenBankProcessor(getabspath(path), config=config, timeout=4)
#     except:
#         logger.exception("Error setting up run_step")
#         return HttpResponse('ERROR', status=500)

#     try:
#         pspath = gbp.process()
#         logger.debug("Finished processing.")
#         pspath = getrelpath(pspath)
#         #url = reverse('results', args=[psname]) + encode_config(config, path=path)
#         result = {'next':'results',
#                   'download_url': get_raw_url(request, pspath),
#                   'steps': gbp.timer.steps}

#     except Timeout, pt:
#         result = {'next':'process',
#                   'steps': pt.steps,
#                   }
#     except:
#         logger.exception("Error in run_step")
#         result = {'next':'ERROR', 'steps': gbp.timer.steps}
#         status=500

#     try:
#         ago = config.get('acgt_gamma_output')
#         if ago:
#             result['files'] = [get_raw_url(request, os.path.join(ago,v))
#                                for v in os.listdir(ago)]
#         # files = set([get_raw_url(request, v)
#         #              for (k,v) in config.items()
#         #              if v and (k in gbp.AP_file_keys)])

#     except:
#         logger.exception("Error building file list.")

#     return HttpResponse(json.dumps(result), status=status)


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


def kickstart(request, path, config):

    #TODO: Some logic here to see if they're already running a job and
    #propose cancelling it

    jobid = client.enqueue(main.process_all, [getabspath(path), config])
    return jobid

def send_email(config, results_link, email_address):
    #config is the result of running the process, needs to be first parameter
    try:
        logger.debug("Task completed; sending email to %r", email_address)
        from django.core.mail import EmailMultiAlternatives
        subject = 'NPACT results ready for "{0}"'.format(config['first_page_title'])
        plaintext = get_template('email-results.txt')
        htmly     = get_template('email-results.html')

        d = Context({ 'keep_days': settings.ATIME_DEFAULT,
                      'results_link': results_link})

        text_content = plaintext.render(d)
        html_content = htmly.render(d)
        msg = EmailMultiAlternatives(subject, text_content, to=[email_address])
        msg.attach_alternative(html_content, "text/html")
        msg.send(fail_silently=False)
        logger.debug("Finished sending email.")
    except:
        logger.exception("Failed sending email to %r", to_address)
