from django.conf.urls.defaults import patterns, include, url


urlpatterns = patterns('spat.views',
     url(r'^run', 'run'),
     url(r'^start', 'start.view'),
     url(r'^library', 'library'),
     url(r'^$', 'index'))
