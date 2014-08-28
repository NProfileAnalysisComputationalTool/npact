import os.path
import glob
import logging

from django.conf import settings
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect
from django.template import RequestContext
from django.views import static
from django.contrib import messages

from pynpact import prepare
from npactweb import library_root


logger = logging.getLogger(__name__)


#### Views that are small enough to be inline here.

def static_serve_wrapper(request, path):
    if settings.DEBUG:
        return static.serve(
            request, path=path, document_root=settings.MEDIA_ROOT)


# def library(request):
#     """List all the files in the library (downloaded from NCBI)

#     This is currently unused.
#     """
#     files = []
#     for gbk in glob.iglob(os.path.join(library_root(), "*.gbk")):
#         data = prepare.try_parse(gbk)

#         if not data.get('url'):
#             data['url'] = reverse('run', args=['library/' + data['basename']])

#         files.append(data)

#     logger.debug("Found %d files in library", len(files))
#     return render_to_response('library.html',
#                               {'files': files},
#                               context_instance=RequestContext(request))


def view_none(request):
    messages.error(
        request,
        "No genome source selected, please upload one, "
        "or go to the library and select one.")
    return HttpResponseRedirect(reverse('start'))
