from settings_shared import *

TEMPLATE_DIRS = (
    "/var/www/fooproj/fooproj/templates",
)

MEDIA_ROOT = '/var/www/fooproj/uploads/'
# put any static media here to override app served static media
STATICMEDIA_MOUNTS = (
    ('/sitemedia', '/var/www/fooproj/fooproj/sitemedia'),	
)


DEBUG = False
TEMPLATE_DEBUG = DEBUG

try:
    from local_settings import *
except ImportError:
    pass
