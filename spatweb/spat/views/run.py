# Create your views here.
import os.path, logging, tempfile

from django.conf import settings
from django.http import HttpResponse,HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.template import RequestContext
from django import forms
from django.contrib import messages

from ordereddict import OrderedDict

from __init__ import session_key


from pynpact import prepare

# Get an instance of a logger
logger = logging.getLogger(__name__)


class RunForm(forms.Form) :
    first_page_title = forms.CharField()
    following_page_title = forms.CharField(required=False)


def try_parse(request) :
    data = OrderedDict()
    path = request.session[session_key("input_path")]
    gb_record = None


    try :
        gbrec = prepare.open_parse(path)
        data['length'] = len(gbrec)
        data['id'] = gbr.id
        data['date'] = gbr.annotations.get('date')
        data['description'] = gbr.description
        return data

    except :
        return None

def prefill_form(request) :
    if request.method == 'POST' :
        return RunForm(request.POST)
    else :
        data = try_parse(request)
        form = RunForm()
        if data :
            form.first_page_title.initial = data.get['description']
            form.following_page_title.initial = data.get['description']


def view(request):
    if not session_key("input_path") in request.session :
        messages.error(request, "No genome source selected, please upload one, or go to the library and select one.")
        return HttpResponseRedirect(reverse('spat.views.start.view'))


    form = prefill_form(request)
    return render_to_response('run.html',{},
                               context_instance=RequestContext(request))
