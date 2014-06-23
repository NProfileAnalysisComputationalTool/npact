import logging

from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.core.exceptions import PermissionDenied
from django.http import HttpResponseRedirect
from django.template import RequestContext

from npactweb import management
from taskqueue import client, tqdaemon


# Get an instance of a logger
logger = logging.getLogger(__name__)


def view(request):
    if not(request.user.is_authenticated() or settings.DEBUG):
        raise PermissionDenied()

    if request.method == 'POST':
        handle_post(request)
        # want to redirect so that refresh works again.
        return HttpResponseRedirect(reverse('management'))
    daemon_status = 'running' if tqdaemon.status() else 'stopped'
    return render_to_response('management.html',
                              {'settings': settings},
                              context_instance=RequestContext(request))


def handle_post(request):
    action = request.POST.get('action')
    if action == 'start-daemon':
        start_daemon(request)
    elif action == 'restart-daemon':
        stop_daemon(request)
        start_daemon(request)
    elif action == 'kill-daemon':
        count = tqdaemon.kill()
        messages.info(request, "Killed {0} processes".format(count))
    elif action == 'cleanup':
        cleanup(request)
    elif action == 'clear-library':
        clear_library(request)

    messages.info(request, "Handled {0}".format(action))


def start_daemon(request):
    client.ensure_daemon()
    if tqdaemon.status():
        messages.info(request, "Daemon started successfully.")
    else:
        messages.error(request, "Daemon failed to start.")


def stop_daemon(request):
    rc = tqdaemon.stop()
    if rc == 0:
        messages.info(request, "Daemon stopped successfully.")


def cleanup(request):
    try:
        days = int(request.POST.get('days', settings.ATIME_DEFAULT))
        if management.cleanup_old_files(days):
            messages.info(request,
                          "Successfully purged files older than %d days." % days)
        else:
            messages.error(request, "Error removing old files.")
    except Exception as e:
        messages.error(request, "Error removing old files: %s" % e)

    try:
        stdout, stderr = management.report_file_size()
        if stdout:
            for l in stdout.split('\n'):
                messages.info(request, l)
        if stderr:
            messages.error(request, stderr)
    except Exception as e:
        messages.error(request, 'Error finding file size' % e)


def clear_library(request):
    try:
        management.clear_library()
    except Exception as e:
        messages.error(request, "Error clearing library: %s" % e)
