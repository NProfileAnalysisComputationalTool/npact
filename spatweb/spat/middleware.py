from django.http import HttpResponse,HttpResponseRedirect

class AddItemsDict(object) :
    def process_request(self, request) :
        request.items = {}
        return None

class RedirectException(Exception):
    def __init__(self, url):
        self.url = url


class RedirectExceptionHandler(object) :
    def process_exception(self, request, exception):
        if isinstance(exception, RedirectException):
            return HttpResponseRedirect(exception.url)
