from django.conf.urls.defaults import patterns, include, url
from django.conf import settings


spatpatterns = \
            patterns('spatweb.views',
                     url(r'^$', 'index', name='index'),
                     url(r'^start$', 'start.view', name="start"),
                     url(r'^start/efetch/(\d+)', 'start.efetch', name="efetch"),
                     url(r'^library', 'library'),
                     url(r'^run/(.+\.gbk?)', 'run.view', name="run"),
                     url(r'^run^', 'run.view_none'),
                     url(r'^results/(.+)', 'run.results', name='results'),
                     )


if settings.DEBUG :
    spatpatterns += patterns('', url(r'^raw/(?P<path>.*)$', 'django.views.static.serve', {
        'document_root': settings.MEDIA_ROOT,
        }, name="raw"))


urlpatterns = patterns('', ('^spat/',include(spatpatterns)))
