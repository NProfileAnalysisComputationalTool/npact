import os.path, glob, logging

from ordereddict import OrderedDict

from django.conf import settings
from django.core.urlresolvers import reverse
from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.template import RequestContext
from django import forms
from django.contrib import messages

from django.views.static import serve

from pynpact import prepare

logger = logging.getLogger(__name__)

#### Helper functions used by all the views

#the session_key for all of spat's data.
def session_key(part) :
    return "spat." + part

def library_root() :
    return os.path.join(settings.MEDIA_ROOT, 'library')


def is_clean_path(path) :
    return path.find('..') < 0


##### Views that are small enough to be inline here.

def index(request) :
    return render_to_response('index.html',{},
                               context_instance=RequestContext(request))

def library(request) :
    files = []
    for gbk in glob.iglob(os.path.join(library_root(),"*.gbk")) :
        data = prepare.try_parse(gbk)

        if not data.get('url') :
            data['url'] =  reverse('run', args=['library/' + data['basename']])

        files.append(data)

    logger.debug("Found %d files in library", len(files))
    return render_to_response('library.html',
                              {'files':files},
                              context_instance=RequestContext(request))


