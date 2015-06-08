from django.conf.urls.defaults import patterns, include, url
from django.conf import settings
from django.views.generic.simple import direct_to_template


npact_patterns = \
    patterns('npactweb.views',
             url(r'^$', 'start.view', name='start'),
             url(r'^about$', direct_to_template,
                 {'template': 'about.html'}, name='about'),
             url(r'^downloads$', direct_to_template,
                 {'template': 'downloads.html'}, name='downloads'),
             url(r'^efetch/(\d+)', 'start.efetch', name="efetch"),
             url(r'^run/(.+)', 'run.run_frame', name="run"),
             url(r'^runstatus/(.*)', 'run.run_status', name='runstatus'),
             url(r'^kickstart/(.*)', 'run.kickstart', name='kickstart'),
             url(r'^translate', 'run.translate', name='translate'),
             url(r'^acgt_gamma_file_list/(.*)', 'run.acgt_gamma_file_list', name='acgt_gamma_file_list'),
             url(r'^(run|config)^', 'view_none'),
             url(r'^raw/(?P<path>.*)$', 'static_serve_wrapper', name='raw'),
             url(r'^management', 'management.view', name='management'))

urlpatterns = patterns('', ('^npact/', include(npact_patterns)))

from django.contrib.staticfiles.urls import staticfiles_urlpatterns
urlpatterns += staticfiles_urlpatterns()
