from django.conf.urls.defaults import patterns, include, url
from django.conf import settings
from django.views.generic.simple import direct_to_template


spatpatterns = \
            patterns('spatweb.views',
                     url(r'^$' , 'start.view', name='start'),
                     url(r'^about$' , direct_to_template, {'template': 'about.html'}, name='about'),
                     url(r'^efetch/(\d+)' , 'start.efetch', name="efetch"),
                     url(r'^library' , 'library'),
                     url(r'^config/(.+\.gbk?)', 'run.config', name="config"),
                     url(r'^run/(.+\.gbk?)' , 'run.run', name="run"),
                     url(r'^(run|config)^' , 'view_none'),
                     url(r'^results/(.+)' , 'run.results', name='results'),
                     url(r'^raw/(?P<path>.*)$', 'static_serve_wrapper', name='raw'),
                     )

urlpatterns = patterns('', ('^spat/', include(spatpatterns)))
