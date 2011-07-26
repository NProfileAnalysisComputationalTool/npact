# Create your views here.
import os.path, logging, tempfile

from django.conf import settings
from django.http import HttpResponse,HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.template import RequestContext
from django import forms
from django.contrib import messages

from __init__ import session_key

from pynpact import main

# Get an instance of a logger
logger = logging.getLogger(__name__)

def view(request):
    return render_to_response('run.html',{},
                               context_instance=RequestContext(request))
