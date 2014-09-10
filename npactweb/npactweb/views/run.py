# Create your views here.
import logging
import os
import os.path
import json

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

from npactweb import assert_clean_path, getabspath, getrelpath
from npactweb.middleware import RedirectException

from taskqueue import client, NoSuchTaskError

# Get an instance of a logger
logger = logging.getLogger(__name__)


def get_raw_url(request, path):
    # return request.build_absolute_uri(reverse('raw', path))
    # TODO: DEPRECATE in favor of get_result_link
    if path.startswith('/'):
        path = getrelpath(path)
    return reverse('raw', args=[path])


def get_result_link(path):
    if path.startswith('/'):
        path = getrelpath(path)
    return reverse('raw', args=[path])


def get_ti(size):
    return forms.TextInput(attrs={'size': size})


class ConfigForm(forms.Form):
    first_page_title = forms.CharField(widget=get_ti(40))
    following_page_title = forms.CharField(required=False, widget=get_ti(40))
    length = forms.IntegerField(required=True, min_value=0,
                                widget=get_ti(8))
    nucleotides = forms.MultipleChoiceField(
        choices=[(i, i) for i in ['a', 'c', 'g', 't']],
        widget=forms.CheckboxSelectMultiple())
    skip_prediction = forms.BooleanField(required=False)
    significance = forms.ChoiceField(
        choices=prepare.significance_levels, required=False,
        label="Prediction Significance")
    start_base = forms.IntegerField()
    end_base = forms.IntegerField()
    alternate_colors = forms.BooleanField(
        required=False, label="Alternative colors")
    email = forms.EmailField(required=False)

    def clean(self):
        cleaned_data = self.cleaned_data
        start = cleaned_data.get('start_base')
        end = cleaned_data.get('end_base')
        if start and end and start > end:
            raise forms.ValidationError(
                "End page must be greater than or equal to start page.")
        return cleaned_data


def get_display_items(request, config):
    yield ('Filename', config['basename'])
    for key in ['date', 'length', 'description']:
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
            url = reverse('run', args=[path]) \
                  + "?" + urlencode(form.cleaned_data, True)
            return HttpResponseRedirect(url)
    else:
        form = ConfigForm(initial=config)

    for key, field in form.fields.items():
        if key in prepare.CONFIG_HELP_TEXT:
            field.help_text = prepare.CONFIG_HELP_TEXT[key]
        elif settings.DEBUG:
            logger.error("Help text missing for config form field: %r", key)

    return render_to_response(
        'config.html', {'form': form, 'parse_data': config,
                        'def_list_items': get_display_items(request, config)},
        context_instance=RequestContext(request))


# Variables, mostly for debugging, that aren't exposed anywhere but
# that if present in the request should be added to the config
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
        logger.exception(
            "Error parsing gbk: %r", getabspath(path, raise_on_missing=False))
        messages.error(request,
                       "There was a problem loading file '%s', "
                       "please try again or try a different record." % path)
        raise RedirectException(reverse('start'))

    cf = ConfigForm(request.REQUEST)
    for f in cf.visible_fields():
        try:
            v = f.field.clean(f.field.to_python(f.data))
            if v:
                logger.debug("Including %r:%r from request.", f.name, v)
                config[f.name] = v
        except forms.ValidationError:
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
    """This is the main processing page.

    Rather this is the frame of the main processing page, a lot of the
    detail is done in javascript on the client.

    In here we kickstart the computation and serve down set of jobids
    for the client to work with.

    """
    assert_clean_path(path, request)

    return render_to_response(
        'processing.html',
        {'status_base': reverse('runstatus', args=['']),
         'kickstart_base': reverse('kickstart', args=[path]),
         'fetch_base': reverse('raw', args=['']),
         'acgt_gamma_base': reverse('acgt_gamma_file_list', args=[''])},
        context_instance=RequestContext(request))


def kickstart(request, path):
    assert_clean_path(path, request)
    config = build_config(path, request)

    # at the moment we always want PDF
    config['pdf_output'] = True

    reconfigure_url = reverse('config', args=[path])
    reconfigure_url += "?" + urlencode_config(config)
    config['reconfigure_url'] = reconfigure_url
    email = request.GET.get('email')
    client.ensure_daemon()
    config = main.process('allplots', config, executor=client.get_server())
    if email:
        target_file = results.pdf_filename or results.combined_ps_name
        assert target_file, \
            "Configured for email but didn't get emailable file."
        # The direcect download link for the PS or PDF file.
        result_link = request.build_absolute_uri(
            get_result_link(target_file))
        # build path back to run screen.
        run_link = reverse('run', args=[path])
        run_link += '?' + urlencode_config(config, exclude=['email'])
        run_link = request.build_absolute_uri(run_link)

        task = util.Task(send_email, email, config, run_link, result_link)
        eid = client.enqueue(task, after=[target_file])
        logger.info("Scheduled email to %r with jobid: %s", email, eid)

    return HttpResponse(
        json.dumps(sanitize_config_for_client(config)),
        content_type="application/json")


def sanitize_config_for_client(config):
    output = {}
    for k, v in config.iteritems():
        if isinstance(v, basestring) and v.startswith(settings.MEDIA_ROOT):
            v = v[len(settings.MEDIA_ROOT):]
            v = v.lstrip("/")
        output[k] = v
    del output['psnames']
    return output


def send_email(email_address, config, run_link, result_link):
    # config is the result of running the process, needs to be first parameter
    try:
        logger.debug("Task completed; sending email to %r", email_address)
        from django.core.mail import EmailMultiAlternatives
        subject = 'NPACT results ready for "{0}"'.format(
            config['first_page_title'])
        plaintext = get_template('email-results.txt')
        htmly = get_template('email-results.html')

        d = Context({'keep_days': settings.ATIME_DEFAULT,
                     'result_link': result_link,
                     'run_link': run_link})

        text_content = plaintext.render(d)
        html_content = htmly.render(d)
        msg = EmailMultiAlternatives(subject, text_content, to=[email_address])
        msg.attach_alternative(html_content, "text/html")
        msg.send(fail_silently=False)
        logger.debug("Finished sending email.")
    except:
        logger.exception("Failed sending email to %r", email_address)


def run_status(request, jobid):
    "This checks on the status of jobid"
    result = {'tid': jobid}
    status = 200
    try:
        result['ready'] = client.ready(getabspath(jobid, raise_on_missing=False))
    except NoSuchTaskError:
        result['message'] = 'Unknown task identifier.'
        status = 404
    except Exception:
        result['message'] = 'Fatal error.'
        status = 500
        logger.exception("Error getting job status. %r", jobid)

    return HttpResponse(json.dumps(result),
                        status=status, content_type="application/json")


def acgt_gamma_file_list(request, acgt_gamma_output):
    acgt_gamma_output = getabspath(acgt_gamma_output)
    files = map(getrelpath, acgt_gamma_output.listdir())
    return HttpResponse(json.dumps(files),
                        status=200, content_type="application/json")
