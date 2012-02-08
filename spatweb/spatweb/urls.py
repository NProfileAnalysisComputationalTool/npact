from django.conf.urls.defaults import patterns, include, url
from django.conf import settings
from django.views.generic.simple import direct_to_template


spatpatterns = \
            patterns('spatweb.views',
                     url(r'^$' , 'start.view', name='start'),
                     url(r'^about$' , direct_to_template, {'template': 'about.html'}, name='about'),
                     url(r'^start/efetch/(\d+)' , 'start.efetch', name="efetch"),
                     url(r'^library' , 'library'),
                     url(r'^run/(.+\.gbk?)' , 'run.view', name="run"),
                     url(r'^run^' , 'run.view_none'),
                     url(r'^results/(.+)' , 'run.results', name='results'),
                     )


if settings.DEBUG :
    spatpatterns += patterns('', url(r'^raw/(?P<path>.*)$', 'django.views.static.serve', {
        'document_root': settings.MEDIA_ROOT,
        }, name="raw"))


urlpatterns = patterns('', ('^spat/', include(spatpatterns)))
