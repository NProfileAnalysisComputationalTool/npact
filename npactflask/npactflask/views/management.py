from flask import url_for, redirect, request, flash, render_template
from werkzeug.exceptions import Unauthorized
from npactflask import management, app
from taskqueue import client, tqdaemon


# Get an instance of a logger
logger = app.logger


def view():
    if (not app.config['DEBUG'] and request.environ.get('REMOTE_USER') is None):
        raise Unauthorized()
    if request.method == 'POST':
        handle_post()
        # want to redirect so that refresh works again.
        return redirect(url_for('management'))
    daemon_status = 'running' if tqdaemon.status() else 'stopped'
    return render_template('management.html',
                           **{'app': app,
                              'daemon_status': daemon_status})


def handle_post():
    action = request.form.get('action')
    if action == 'start-daemon':
        start_daemon()
    elif action == 'restart-daemon':
        stop_daemon()
        start_daemon()
    elif action == 'kill-daemon':
        count = tqdaemon.kill()
        flash(request, "Killed {0} processes".format(count))
    elif action == 'cleanup':
        cleanup()
    elif action == 'clear-library':
        clear_library()

    flash(request, "Handled {0}".format(action))


def start_daemon():
    client.ensure_daemon()
    if tqdaemon.status():
        flash("Daemon started successfully.")
    else:
        flash("Daemon failed to start.")


def stop_daemon():
    rc = tqdaemon.stop()
    if rc == 0:
        flash("Daemon stopped successfully.")


def cleanup():
    try:
        days = int(request.form.get('days', app.config['ATIME_DEFAULT']))
        if management.cleanup_old_files(days):
            flash("Successfully purged files older than %d days." % days)
        else:
            flash("Error removing old files.")
    except Exception as e:
        flash("Error removing old files: %s" % e)

    try:
        stdout, stderr = management.report_file_size()
        if stdout:
            for l in stdout.split('\n'):
                flash(l)
        if stderr:
            flash(stderr)
    except Exception as e:
        flash('Error finding file size' % e)


def clear_library():
    try:
        management.clear_library()
    except Exception as e:
        flash("Error clearing library: %s" % e)
