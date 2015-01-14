# Create your views here.
import logging
import os
import os.path
import json

from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from django.utils.http import urlencode
from django.template import RequestContext, Context
from django.template.loader import get_template

from pynpact import parsing, util
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

# Variables, mostly for debugging, that aren't exposed anywhere but
# that if present in the request should be added to the config
MAGIC_PARAMS = ['raiseerror', 'force']

VALID_KEYS = ('first_page_title', 'following_page_title', 'nucleotides',
              'significance', 'alternate_colors', 'startBase', 'endBase')


def build_config(path, request):
    "Tries to build the config dictionary for the given path"
    assert_clean_path(path, request)

    try:
        config = parsing.initial(getabspath(path))
    except:
        logger.exception(
            "Error parsing gbk: %r", getabspath(path, raise_on_missing=False))
        messages.error(request,
                       "There was a problem loading file '%s', "
                       "please try again or try a different record." % path)
        raise RedirectException(reverse('start'))

    parsing.detect_format(config)
    for k in VALID_KEYS:
        if k in request.GET and request.GET[k]:
            config[k] = request.GET[k]
    # These need to be converted to ints
    parsing.startBase(config)
    parsing.endBase(config)
    parsing.first_page_title(config)



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
    keys = VALID_KEYS + MAGIC_PARAMS
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
        {
            'status_base': reverse('runstatus', args=['']),
            'kickstart_base': reverse('kickstart', args=[path]),
            'fetch_base': reverse('raw', args=['']),
            'acgt_gamma_base': reverse('acgt_gamma_file_list', args=['']),
            'base_href': reverse('run', args=[path])
        },
        context_instance=RequestContext(request))


def kickstart(request, path):
    assert_clean_path(path, request)
    config = build_config(path, request)
    verb = request.GET['verb']
    verb = verb.split(',')

    email = request.GET.get('email')
    client.ensure_daemon()
    executor = client.get_server()

    for v in verb:
        if email:
            config = main.process('allplots', config, executor=executor)
            build_email(request, path, config)
        elif v == 'parse':
            # we've already parsed in `build_config` above.
            pass
        else:
            # main.process handles the rest of the verbs, or errors
            config = main.process(v, config, executor=executor)

    return HttpResponse(
        json.dumps(sanitize_config_for_client(config)),
        content_type="application/json")


def build_email(request, path, config):
    try:
        email = request.GET.get('email')
        target_file = config['pdf_filename'] or config['combined_ps_name']
        assert target_file, \
            "Configured for email but didn't get emailable file."
        # The direct download link for the PS or PDF file.
        result_link = request.build_absolute_uri(
            get_result_link(target_file))
        # build path back to run screen.
        run_link = reverse('run', args=[path])
        run_link += '?' + urlencode_config(config, exclude=['email'])
        run_link = request.build_absolute_uri(run_link)

        task = util.Task(send_email, email, config, run_link, result_link)
        eid = client.enqueue(task, after=[target_file])
        logger.info("Scheduled email to %r with jobid: %s", email, eid)
    except:
        logger.exception("Error scheduling email to send")


def sanitize_config_for_client(config):
    output = {}
    for k, v in config.iteritems():
        if isinstance(v, basestring) and v.startswith(settings.MEDIA_ROOT):
            v = v[len(settings.MEDIA_ROOT):]
            v = v.lstrip("/")
        output[k] = v
    if 'psnames' in output:
        del output['psnames']
    return output


def send_email(email_address, config, run_link, result_link):
    try:
        logger.debug("Task completed; sending email %r, %r, %r",
                     email_address, run_link, result_link)
        from django.core.mail import EmailMultiAlternatives
        subject = 'NPACT results ready for "{0}"'.format(
            config['first_page_title'])
        plaintext = get_template('email-results.txt')
        htmly = get_template('email-results.html')

        d = Context({'keep_days': settings.ATIME_DEFAULT,
                     'results_link': result_link,
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
        abspath = getabspath(jobid, raise_on_missing=False)
        # Task IDs == filenames on disk, if the file exists don't even
        # bother to check with the server
        if os.path.exists(abspath):
            result['ready'] = True
        elif client.ready(abspath):
            result['ready'] = True
            # Ensure the task finished correctly
            client.result(abspath)
        else:
            result['ready'] = False

    except NoSuchTaskError:
        result['message'] = 'Unknown task identifier.'
        status = 404
    except Exception as e:
        result['message'] = 'Fatal error.'
        result['exception'] = str(e)
        status = 500
        logger.exception("Error getting job status. %r", jobid)

    return HttpResponse(json.dumps(result),
                        status=status, content_type="application/json")


def acgt_gamma_file_list(request, acgt_gamma_output):
    acgt_gamma_output = getabspath(acgt_gamma_output)
    files = map(getrelpath, acgt_gamma_output.listdir())
    return HttpResponse(json.dumps(files),
                        status=200, content_type="application/json")
