import logging
import os
import os.path
import json
from datetime import datetime

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

from taskqueue import client, NoSuchTaskError


# Get an instance of a logger
logger = logging.getLogger(__name__)

def view(request):
    if not(request.user.is_authenticated() or settings.DEBUG):
        raise PermissionDenied()
        
    return render_to_response('management.html',
                              {},
                              context_instance=RequestContext(request))
