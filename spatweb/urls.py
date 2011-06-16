from django.conf.urls.defaults import patterns, include, url
from django.conf.urls.defaults import handler404

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
     url(r'.', 'spat.views.run'),
     url(r'^$', 'spat.views.index'),
    # url(r'^spatweb/', include('spatweb.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # url(r'^admin/', include(admin.site.urls)),
)
