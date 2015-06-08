import logging
import sys
from datetime import datetime

from django.conf import settings
from django.views.debug import technical_500_response
from django.http import HttpResponseRedirect, HttpResponse
from npactweb import MissingFileError, ImmediateHttpResponseException, \
    RedirectException


logger = logging.getLogger(__name__)


class AddItemsDict(object):
    def process_request(self, request):
        request.items = {}
        return None


class NPactResponseExceptionHandler(object):
    def process_exception(self, request, exception):
        if isinstance(exception, ImmediateHttpResponseException):
            return exception.httpResponse
        elif isinstance(exception, RedirectException):
            return HttpResponseRedirect(exception.url)
        elif isinstance(exception, MissingFileError):
            return HttpResponse(exception.message, status=404,
                                content_type="text/plain")
        else:
            logger.exception(exception)
            return HttpResponse(repr(exception), status=500,
                                content_type="application/json")


class SaveTraceExceptionMiddleware(object):
    """
     - add to middleware 'SaveTraceExceptionMiddleware',
    """
    def process_exception(self, request, exception):
        if not settings.DEBUG and settings.EXC_TRAC_PATH:
            try:
                tech_response = technical_500_response(
                    request, *sys.exc_info())

                error_id = datetime.now().isoformat("-").replace(":", "-")
                fname = "{0}/{1}.txt".format(settings.EXC_TRACE_PATH, error_id)
                with open(fname, "w") as fout:
                    fout.write(tech_response.content)
                logger.info(
                    "Exception technical response saved in %s", fout.name)
            except Exception, e:
                logger.error(
                    "Error when saving exception to file: '%s' / '%s' ",
                    str(sys.exc_info()), str(e))
