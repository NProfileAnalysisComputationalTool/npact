# Create your views here.
from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.template import RequestContext

from django.contrib import messages


def index(request) :
    return render_to_response('index.html',{},
                               context_instance=RequestContext(request))

def library(request) :
    return render_to_response('library.html', {},
                               context_instance=RequestContext(request))

def start(request) :
    return render_to_response('start.html',{},
                               context_instance=RequestContext(request))

def run(request):
    messages.info(request,"Method was %s" % request.method)
    return render_to_response('run.html',{},
                               context_instance=RequestContext(request))

def results(request) :
    pass

