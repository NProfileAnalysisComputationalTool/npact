import os.path, glob, logging

from django.conf import settings
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.utils.http import urlencode
from django.views import static

from pynpact import prepare
from npactweb import library_root


logger = logging.getLogger(__name__)


##### Views that are small enough to be inline here.

def static_serve_wrapper(request, path):
    if settings.DEBUG:
        return static.serve(request, path=path, document_root=settings.MEDIA_ROOT)
        

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


def view_none(request) :
    messages.error(request,
                   "No genome source selected, please upload one, "
                   "or go to the library and select one.")
    return HttpResponseRedirect(reverse('start'))
