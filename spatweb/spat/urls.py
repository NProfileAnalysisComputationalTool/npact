from django.conf.urls.defaults import patterns, include, url


urlpatterns = patterns('spat.views',
     url(r'^run/(.+\.gbk?)', 'run.view', name="run"),
     url(r'^run^', 'run.view_none'),
     url(r'^start', 'start.view', name="start"),
     url(r'^library', 'library'),
     url(r'^$', 'index'))
