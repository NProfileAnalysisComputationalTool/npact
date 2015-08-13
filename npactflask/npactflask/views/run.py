import flask
import os
import os.path
import Bio.Seq

from path import path as Path
from flask import (
    url_for, request, flash, redirect, json, jsonify,
    send_from_directory, send_file
)
from werkzeug.exceptions import NotFound
from pynpact import main, parsing, util, executors
from npactflask import app
from npactflask.views import getabspath, getrelpath

gexec = executors.GeventExecutor()


logger = app.logger
VALID_KEYS = ('first_page_title', 'following_page_title', 'nucleotides',
              'significance', 'alternate_colors', 'startBase', 'endBase',
              'basesPerGraph', 'x-tics', 'mycoplasma')
MAGIC_PARAMS = ('raiseerror', 'force')


def get_result_link(path):
    return url_for('raw', path=path)


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
        config['nucleotides'] = request.args.getlist('nucleotides')

    for key in MAGIC_PARAMS:
        if key in request.args:
            config[key] = request.args[key]

    return config


def dicforurl(config, exclude=None):
    """Reduce a config dictionary suitable for passing in the url
    """
    keys = VALID_KEYS + MAGIC_PARAMS
    if exclude:
        keys = set(keys) - set(exclude)
    return util.reducedict(config, keys)


def run_frame(path):
    """This is the main processing page.

    Rather this is the frame of the main processing page, a lot of the
    detail is done in javascript on the client.

    In here we kickstart the computation and serve down set of jobids
    for the client to work with.

    """

    if not os.path.exists(getabspath(path, False)):
        flash("Couldn't find genome in: %s" % path)
        return redirect(url_for('start'))

    return flask.render_template(
        'processing.html', **{
            'BASE_URL': app.config['APPLICATION_ROOT'],
            'PATH': path,
            'status_base': url_for('.run_status', path=''),
            'kickstart': url_for('kickstart', path=path),  # args=[path]
            'translate': url_for('translate'),
            'fetch_base': url_for('raw', path=''),
            'acgt_gamma_base': url_for('.acgt_gamma_file_list', path=''),
            'STATIC_BASE': url_for('static', filename=''),
            'base_href': url_for('run_frame', path=path)
        })


def translate():
    try:
        # table 4 is for mycoplasma ala:
        # http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
        table = 1
        if request.form.get('mycoplasma'):
            table = 4
        seq = Bio.Seq.Seq(request.form.get('seq'))
        rc = int(request.form.get('complement', 0))
        if rc:
            seq = seq.reverse_complement()
        trans = Bio.Seq.translate(seq, table)
        return jsonify({
            'seq': str(trans),
            'complement': rc and str(seq)})
    except Exception as e:
        logger.exception('Error While Translating')
        return flask.make_response(repr(e), 500)


def kickstart(path):
    config = build_config(path)
    verb = request.args['verb']
    verb = verb.split(',')
    email = request.args.get('email')
    for v in verb:
        if email:
            config = main.process('allplots', config, executor=gexec)
            build_email(path, config)
        elif v == 'parse':
            # we've already parsed in `build_config` above.
            pass
        else:
            # main.process handles the rest of the verbs, or errors
            config = main.process(v, config, executor=gexec)

    return flask.make_response(
        json.dumps(sanitize_config_for_client(config)))


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
        run_link = url_for(
            'run', path=path,
            **dicforurl(config, exclude=['email']))
        run_link = request.build_absolute_uri(run_link)

        task = util.Task(send_email, email, config, run_link, result_link)
        eid = gexec.enqueue(task, after=[target_file])
        logger.info("Scheduled email to %r with jobid: %s", email, eid)
    except:
        logger.exception("Error scheduling email to send")


def sanitize_config_for_client(config):
    output = {}
    uploads = app.config['UPLOADS']
    for k, v in config.iteritems():
        if isinstance(v, basestring) and v.startswith(uploads):
            v = v[len(uploads):]
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

        d = app.app_context({'keep_days': app.config['ATIME_DEFAULT'],
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
        elif gexec.ready(abspath):
            result['ready'] = True
            # Ensure the task finished correctly, this could except
            gexec.result(abspath)
        else:
            result['ready'] = False

    except KeyError:  # TODO: Use a restored NoSuchTak
        result['message'] = 'Unknown task identifier.'
        status = 404
    except Exception as e:
        result['message'] = 'Fatal error.'
        result['exception'] = str(e)
        status = 500
        logger.exception("Error getting job status. %r", path)

    return flask.make_response(json.dumps(result), status)


@app.route('/blockon/<path:path>')
def blockon(path):
    # Task IDs == filenames on disk
    abspath = getabspath(path, raise_on_missing=False)
    try:
        abspath = getabspath(path, raise_on_missing=False)
        # Task IDs == filenames on disk, if the file exists don't even
        # bother to check with the server
        if os.path.exists(abspath) or gexec.result(abspath, timeout=None):
            return send_file(abspath)
        else:
            raise NotFound()
    except KeyError:
        raise NotFound()
    except Exception as e:
        logger.exception("Error getting job status. %r", path)
        r = jsonify(message='Fatal error.', exception=str(e))
        r.status_code = 500
        return r


@app.route('/getpdf/<path:path>')
def getpdf(path):
    config = build_config(path)
    config = main.process('allplots', config, executor=gexec)
    pdf = config['pdf_filename']
    gexec.result(pdf, timeout=None)
    return send_file(pdf)


def acgt_gamma_file_list(path):
    acgt_gamma_output = getabspath(path)
    files = map(getrelpath, acgt_gamma_output.listdir())
    return flask.make_response(json.dumps(files), 200)


@app.route('/acgt_gamma/<path:path>')
def acgt_gamma(path):
    config = build_config(path)
    config = main.process('acgt_gamma', config, executor=gexec)
    tid_output_directory = config['acgt_gamma_output']
    gexec.result(tid_output_directory, timeout=None)
    files = map(getrelpath, Path(tid_output_directory).listdir())
    return jsonify(config=config, files=files)
