import logging
import os
import os.path
import json
from datetime import datetime
import time

from django import forms
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from django.utils.http import urlencode
from django.template import RequestContext, Context
from django.template.loader import get_template
from django.core.exceptions import PermissionDenied

from django.contrib.auth.decorators import login_required

from npactweb import assert_clean_path, getabspath, getrelpath
from npactweb.middleware import RedirectException

from taskqueue import client, tqdaemon


# Get an instance of a logger
logger = logging.getLogger(__name__)

def view(request):
    if not(request.user.is_authenticated() or settings.DEBUG):
        raise PermissionDenied()

    if request.method == 'POST':
        handle_post(request)
        #want to redirect so that refresh works again.
        return HttpResponseRedirect(reverse('management'))
    daemon_status = 'running' if tqdaemon.status() else 'stopped'
    return render_to_response('management.html',
                              locals(),
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
