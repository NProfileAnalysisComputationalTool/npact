from flask import url_for, redirect, request, flash, render_template
from werkzeug.exceptions import Unauthorized
from npactflask import management, app


# Get an instance of a logger
logger = app.logger


def view():
    if not app.config['DEBUG'] and request.environ.get('REMOTE_USER') is None:
        raise Unauthorized()
    if request.method == 'POST':
        handle_post()
        # want to redirect so that refresh works again.
        return redirect(url_for('management'))
    return render_template('management.html',
                           default_atime=app.config['ATIME_DEFAULT'])


def handle_post():
    action = request.form.get('action')
    logger.info("Handling action %r", action)
    if action == 'cleanup':
        cleanup()
    elif action == 'clear-library':
        clear_library()


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
            for l in stdout.strip().split('\n'):
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
