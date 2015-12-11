import flask
import os
import re
import datetime
import os.path
import Bio.Seq

from path import path as Path
from flask import (
    url_for, request, flash, redirect, json, jsonify, send_file,
    render_template
)
from flask_mail import Message, email_dispatched
from werkzeug.exceptions import NotFound
from pynpact import main, parsing, util, executors
from npactflask import app, mail
from npactflask.views import getabspath, getrelpath

gexec = executors.GeventExecutor()


logger = app.logger
VALID_KEYS = ('first_page_title', 'following_page_title', 'nucleotides',
              'significance', 'alternate_colors', 'startBase', 'endBase',
              'basesPerGraph', 'x-tics', 'mycoplasma')
MAGIC_PARAMS = ('raiseerror', 'force')


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
    parsing.mycoplasma(config)

    # fixup nucleotides list
    if 'nucleotides' in request.args:
        config['nucleotides'] = request.args.getlist('nucleotides')

    for key in MAGIC_PARAMS:
        if key in request.args:
            config[key] = request.args[key]

    return config


def dictforurl(config, exclude=None):
    """Reduce a config dictionary suitable for passing in the url
    """
    keys = VALID_KEYS + MAGIC_PARAMS
    if exclude:
        keys = set(keys) - set(exclude)
    return util.reducedict(config, keys)


@app.route('/run/<path:path>')
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

    return render_template(
        'processing.html', **{
            'BASE_URL': app.config['APPLICATION_ROOT'],
            'PATH': path,
            'status_base': url_for('.run_status', path=''),
            'kickstart': url_for('kickstart', path=path),  # args=[path]
            'translate': url_for('translate'),
            'save_track': url_for('save_track'),
            'fetch_base': url_for('raw', path=''),
            'STATIC_BASE': url_for('static', filename=''),
            'base_href': url_for('run_frame', path=path)
        })


@app.route('/translate', methods=['POST'])
def translate():
    "Translate arbitrary sequence"
    try:
        # table 4 is for mycoplasma ala:
        # http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
        table = 1
        config = request.get_json()
        if parsing.mycoplasma(config):
            app.logger.debug('using mycoplasma table=4')
            table = 4
        seq = Bio.Seq.Seq(config.get('seq'))
        rc = parsing.tobool(config.get('complement'))
        if rc:
            seq = seq.reverse_complement()
        trans = Bio.Seq.translate(seq, table)
        return jsonify({
            'seq': str(seq),
            'trans': str(trans)})
    except Exception as e:
        logger.exception('Error While Translating')
        return flask.make_response(repr(e), 500)


def _track_line_from_orf(orf):
    if orf.get('complement'):
        return "{name} complement({start}..{end})\n".format(**orf)
    else:
        return "{name} {start}..{end}\n".format(**orf)


def _upload_root():
    # TODO: whats the actual application root?
    return os.path.join(os.path.abspath(os.curdir), 'webroot/uploads')


def _new_track_file_name(track):
    fn = track.get('filename')
    dt = datetime.datetime.now()
    # strip previous timestamp
    fn = re.sub(r'\.(\d|_)*$', '', fn)
    path = fn + '.' + re.sub(r':|\.|-', '_', dt.isoformat("_"))
    return (os.path.join(_upload_root(), path), path)


@app.route('/save-track', methods=['POST'])
def save_track():
    try:
        config = request.get_json()
        track = config.get('track')
        orfs = track.get('data')
        if not track or not orfs:
            raise ValueError('Must have a track and orfs to save a track')
        (path, rtn) = _new_track_file_name(track)
        with open(path, "w") as f:
            for orf in orfs:
                f.write(_track_line_from_orf(orf))
        return jsonify({'filename': rtn})
    except Exception as e:
        logger.exception('Error While Saving Track')
        return flask.make_response(repr(e), 500)


@app.route('/translate/<path:path>', methods=['GET'])
def translate2(path):
    "Translate the known file with params"
    config = build_config(path)
    try:
        return jsonify(
            parsing.translate(config, parsing.tobool(request.args.get('rc')))
        )
    except Exception as e:
        logger.exception('Error While Translating')
        return flask.make_response(repr(e), 500)


@app.route('/kickstart/<path:path>')
def kickstart(path):
    config = build_config(path)
    verb = request.args['verb']
    verb = verb.split(',')
    email = request.args.get('email')
    for v in verb:
        if email:
            schedule_email(email, path, config)

        elif v == 'parse':
            # we've already parsed in `build_config` above.
            pass
        else:
            # main.process handles the rest of the verbs, or errors
            config = main.process(v, config, executor=gexec)

    return flask.make_response(
        json.dumps(sanitize_config_for_client(config)))


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


@app.route('/runstatus/<path:path>')
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
    pdf = config['allplots_result']
    gexec.result(pdf, timeout=None)
    return send_file(pdf, as_attachment=True)


@app.route('/acgt_gamma/<path:path>')
def acgt_gamma(path):
    config = build_config(path)
    config = main.process('acgt_gamma', config, executor=gexec)
    tid_output_directory = config['acgt_gamma_output']
    gexec.result(tid_output_directory, timeout=None)
    files = map(getrelpath, Path(tid_output_directory).listdir())
    return jsonify(
        NewOrfsFile=getrelpath(config['NewOrfsFile']),
        ModifiedOrfsFile=getrelpath(config['ModifiedOrfsFile']),
        HitsFile=getrelpath(config['HitsFile']),
        files=files)


def schedule_email(to, path, config):
    config = main.process('allplots', config, executor=gexec)
    try:
        target_file = config.get('allplots_result')
        assert target_file, \
            "Configured for email but didn't get emailable file."
        # The direct download link for the PS or PDF file.
        result_link = url_for('raw',
                              path=getrelpath(target_file), _external=True)
        # build path back to run screen.
        run_link = url_for('run_frame', path=path, _external=True,
                           **dictforurl(config, exclude=['email']))
        task = util.Task(send_email, to, path, config, result_link, run_link)
        eid = gexec.enqueue(task, after=[target_file])
        logger.info("Scheduled email to %r with jobid: %s", to, eid)
    except:
        logger.exception("Error scheduling email to send")


def send_email(to, path, config, result_link, run_link):
    try:
        logger.debug("Task completed; sending email %r, %r, %r",
                     to, run_link, result_link)
        subject = 'NPACT results ready for "{0}"'.format(
            config['first_page_title'])

        with app.app_context():
            ctx = {'keep_days': app.config['ATIME_DEFAULT'],
                   'results_link': result_link,
                   'run_link': run_link}
            text_content = render_template('email-results.txt', **ctx)
            html_content = render_template('email-results.html', **ctx)
            msg = Message(subject, recipients=[to])
            msg.body = text_content
            msg.html = html_content
            mail.send(msg)
    except:
        logger.exception("Failed sending email to %r", to)


@email_dispatched.connect_via(app)
def log_mail_finished(message, app):
    logger.debug("Finished sending email to %r", message.recipients)
