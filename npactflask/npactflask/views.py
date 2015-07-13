import flask
import logging
import os
import os.path
import json
import Bio.Seq
from flask import url_for, request, send_from_directory, flash, redirect
from npactflask import app
from npactflask import settings
from pynpact import main, parsing, util
from taskqueue import client, NoSuchTaskError


logger = logging.getLogger(__name__)
VALID_KEYS = ('first_page_title', 'following_page_title', 'nucleotides',
              'significance', 'alternate_colors', 'startBase', 'endBase',
              'basesPerGraph', 'x-tics', 'mycoplasma')
MAGIC_PARAMS = ('raiseerror', 'force')


def library_root():
    abspath = getabspath("library", False)
    if not os.path.exists(abspath):
        logger.info("Creating library root at %s", abspath)
        os.makedirs(abspath)
    return abspath


def getrelpath(abspath):
    return os.path.relpath(abspath, settings.MEDIA_ROOT)


def is_clean_path(path):
    return (not path.startswith('/')) and os.path.normpath(path) == path


def assert_clean_path(path, request,
                      message="Path contained illegal characters. Please "
                      "select a new GBK file.",
                      destination='start'):
    if not is_clean_path(path):
        flash(message)
        raise RedirectException(url_for(destination))


class MissingFileError(Exception):
    pass


class ImmediateHttpResponseException (Exception):
    httpResponse = None

    def __init__(self, httpResponse):
        self.httpResponse = httpResponse


class RedirectException(Exception):
    url = None

    def __init__(self, url):
        self.url = url


def getabspath(relpath, raise_on_missing=True):
    if not is_clean_path(relpath):
        logger.error("Illegal path submitted", relpath)
        raise MissingFileError("Bad characters")
    abspath = os.path.join(settings.MEDIA_ROOT, relpath)
    if raise_on_missing and not os.path.exists(abspath):
        raise MissingFileError("Path not found: " + relpath)
    return abspath


def get_result_link(path):
    return url_for('raw', path=path)


@app.context_processor
def STATIC():
    return dict(STATIC_URL=url_for('static', filename=''))


def index():
    return flask.render_template('start.html')


def about():
    return flask.render_template('about.html')


def view_none():
    flash(
        "No genome source selected, please upload one, "
        "or go to the library and select one.")
    return redirect(url_for('index', filename=''))


def static_serve_wrapper(path):
    if app.debug:
        return send_from_directory(settings.MEDIA_ROOT, path)


def run_frame(path):
    """This is the main processing page.

    Rather this is the frame of the main processing page, a lot of the
    detail is done in javascript on the client.

    In here we kickstart the computation and serve down set of jobids
    for the client to work with.

    """
    return flask.render_template(
        'processing.html', **{
            'status_base': url_for('.run_status', path=''),
            'kickstart_base': url_for('kickstart', path=path),  # args=[path]
            'translate_base': url_for('translate', path=''),
            'fetch_base': url_for('raw', path=''),
            'acgt_gamma_base': url_for('.acgt_gamma_file_list', path=''),
            'base_href': url_for('run_frame', path=path)
        })


def run_status(path):
    "This checks on the status of jobid=path"
    result = {'tid': path}
    status = 200
    try:
        abspath = getabspath(path, raise_on_missing=False)
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
        logger.exception("Error getting job status. %r", path)

    return flask.make_response(json.dumps(result), status)


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


def dicforurl(config, exclude=None):
    """Reduce a config dictionary suitable for passing in the url
    """
    keys = VALID_KEYS + MAGIC_PARAMS
    if exclude:
        keys = set(keys) - set(exclude)
    return util.reducedict(config, keys)


def build_config(path):
    "Tries to build the config dictionary for the given path"
    try:
        config = parsing.initial(getabspath(path))
    except:
        logger.exception(
            "Error parsing gbk: %r", getabspath(path, raise_on_missing=False))
        raise

    parsing.detect_format(config)
    for k in VALID_KEYS:
        v = request.args.get(k)
        if v:
            nv = parsing.number(v)
            config[k] = v if nv is None else nv
    parsing.endBase(config)

    # fixup nucleotides list
    if 'nucleotides' in request.args:
        logger.info("nucleotides: %r", request.args.getlist('nucleotides'))
        config['nucleotides'] = request.args.getlist('nucleotides')

    for key in MAGIC_PARAMS:
        if key in request.args:
            config[key] = request.args[key]

    return config


def build_email(path, config):
    try:
        email = request.args.get('email')
        target_file = config['pdf_filename'] or config['combined_ps_name']
        assert target_file, \
            "Configured for email but didn't get emailable file."
        # The direct download link for the PS or PDF file.
        result_link = request.build_absolute_uri(
            get_result_link(target_file))
        # build path back to run screen.
        run_link = url_for('run', path=path, **dicforurl(config, exclude=['email']))
        run_link = request.build_absolute_uri(run_link)

        task = util.Task(send_email, email, config, run_link, result_link)
        eid = client.enqueue(task, after=[target_file])
        logger.info("Scheduled email to %r with jobid: %s", email, eid)
    except:
        logger.exception("Error scheduling email to send")


def send_email(email_address, config, run_link, result_link):
    try:
        logger.debug("Task completed; sending email %r, %r, %r",
                     email_address, run_link, result_link)
        from django.core.mail import EmailMultiAlternatives
        subject = 'NPACT results ready for "{0}"'.format(
            config['first_page_title'])

        d = app.app_context({'keep_days': settings.ATIME_DEFAULT,
                             'results_link': result_link,
                             'run_link': run_link})

        text_content = flask.render_template('email-results.text', d)
        html_content = flask.render_template('email-results.html', d)
        msg = EmailMultiAlternatives(subject, text_content, to=[email_address])
        msg.attach_alternative(html_content, "text/html")
        msg.send(fail_silently=False)
        logger.debug("Finished sending email.")
    except:
        logger.exception("Failed sending email to %r", email_address)


def kickstart(path):
    config = build_config(path)
    verb = request.args['verb']
    verb = verb.split(',')

    email = request.args.get('email')
    client.ensure_daemon()
    executor = client.get_server()
    for v in verb:
        if email:
            config = main.process('allplots', config, executor=executor)
            build_email(path, config)
        elif v == 'parse':
            # we've already parsed in `build_config` above.
            pass
        else:
            # main.process handles the rest of the verbs, or errors
            config = main.process(v, config, executor=executor)

    return flask.make_response(
        json.dumps(sanitize_config_for_client(config)))


def translate():
    # table 4 is for mycoplasma ala:
    # http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    table = 1
    if request.POST.get('mycoplasma'):
        table = 4
    seq = Bio.Seq.Seq(request.POST.get('seq'))
    rc = request.POST.get('complement')
    if rc:
        seq = seq.reverse_complement()
    trans = Bio.Seq.translate(seq, table)
    return flask.make_response(json.dumps({
        'seq': str(trans),
        'complement': rc and str(seq)}
    ), 200)


def acgt_gamma_file_list(path):
    acgt_gamma_output = getabspath(path)
    files = map(getrelpath, acgt_gamma_output.listdir())
    return flask.make_response(json.dumps(files), 200)


def efetch():
    return
