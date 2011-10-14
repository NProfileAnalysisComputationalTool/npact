import os.path, glob, logging

from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.utils.http import urlencode
from pynpact import prepare
from spatweb import library_root


logger = logging.getLogger(__name__)

#### Helper functions used by all the views

def get_return_url(request):
    if request.GET.get('path'):
        return reverse('run',args=[request.GET.get('path')]) + "?" + urlencode(request.GET)
    else:
        return None

def get_raw_url(request, path):
    return request.build_absolute_uri(reverse('index') + path)
    
        


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

