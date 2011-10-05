from django.conf.urls.defaults import patterns, include, url
from django.conf import settings


urlpatterns = \
            patterns('spat.views',
                     url(r'^$', 'index'),
                     url(r'^start$', 'start.view', name="start"),
                     url(r'^start/efetch/(\d+)', 'start.efetch', name="efetch"),
                     url(r'^library', 'library'),
                     url(r'^run/(.+\.gbk?)', 'run.view', name="run"),
                     url(r'^run^', 'run.view_none'),
                     url(r'^results/(.+)', 'run.results', name='results'),
                     )


if settings.DEBUG :
    urlpatterns += patterns('', url(r'^raw/(?P<path>.*)$', 'django.views.static.serve', {
        'document_root': settings.MEDIA_ROOT,
        }, name="raw"))

else :
    #TODO: "Need to configure a fallback url for apache to punt to."
    raise Exception("Need to configure a fallback url for apache to punt to.")



