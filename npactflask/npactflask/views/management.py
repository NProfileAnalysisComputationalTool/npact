import os
from flask import url_for, redirect, request, flash, render_template, make_response
from werkzeug.exceptions import Unauthorized
from npactflask import cleanup, app
from npactflask.views import getrelpath


# Get an instance of a logger
logger = app.logger



@app.route('/management', endpoint='management', methods=['POST', 'GET'])
def view():
    if request.method == 'POST':
        handle_post()
        # want to redirect so that refresh works again.
        return redirect(url_for('management'))

    return render_template('management.html',
                           environ=request.environ,
                           default_atime=app.config['ATIME_DEFAULT'])


def handle_post():
    action = request.form.get('action')
    logger.info("Handling action %r", action)
    if action == 'cleanup':
        docleanup()


def docleanup():
    try:
        days = int(request.form.get('days', app.config['ATIME_DEFAULT']))
        if cleanup.cleanup_old_files(days):
            flash("Successfully purged files older than %d days." % days, "success")
        else:
            flash("Error removing old files.", "danger")
    except Exception as e:
        flash("Error removing old files: %s" % e, "danger")

    try:
        stdout, stderr = cleanup.report_file_size()
        if stdout:
            for l in stdout.strip().split('\n'):
                size, abspath = l.split('\t', 1)
                relpath = getrelpath(abspath)
                if relpath == '.': relpath = "uploads"
                flash("%s remaining in %s" % (size, relpath), "info")
        if stderr:
            flash(stderr, "danger")
    except Exception as e:
        logger.exception("Error finding size")
        flash('Error finding file size: ' + str(e), "danger")
