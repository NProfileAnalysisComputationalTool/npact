# Create your views here.
from django.http import HttpResponse
from django.shortcuts import render_to_response

def index(request) :
    return render_to_response('index.html',{})

def run(request) :
    return render_to_response('run.html',{})

def library(request) :
    return render_to_response('library.html', {})
