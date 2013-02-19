DEBUG = True
TEMPLATE_DEBUG = DEBUG

import os
from path import path

webroot=(path(__file__).dirname() / "../webroot").realpath()
if not webroot.exists():
    raise Exception("Couldn't find webroot at %s" % webroot)

def ppath(rel, create=False):
    abspath = (webroot / rel).realpath()
    if abspath.exists():
        return abspath
    elif create:
        os.makedirs(abspath)
        return abspath
    else:
        raise Exception("Path '%s' doesn't exist." % abspath)



ADMINS = (
    ('Nathan Bird', 'nathan@acceleration.net'),
    )

MANAGERS = ADMINS

#https://docs.djangoproject.com/en/dev/topics/email/#smtp-backend
EMAIL_USE_TLS = True
EMAIL_HOST = 'smtp.gmail.com'
EMAIL_PORT = 587
EMAIL_HOST_USER = 'npact1.0@gmail.com'
EMAIL_HOST_PASSWORD = 'sictransit2'



TIME_ZONE = 'America/New_York'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = False

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale
USE_L10N = False

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/home/media/media.lawrence.com/media/"
MEDIA_ROOT = ppath('uploads',True)

TQ_DIR = ppath('taskqueue', True)
import taskqueue
taskqueue.BASE_DIR = TQ_DIR
    

# how many days should we keep uploaded files and products that
# haven't been accessed before we delete them.
MEDIA_RETAIN_FOR=7
ATIME_DEFAULT = 30

######## django-mediagenerator settings

MEDIA_DEV_MODE = DEBUG
DEV_MEDIA_URL = '/npact/devassets/'
PRODUCTION_MEDIA_URL = '/npact/assets/'
#GLOBAL_MEDIA_DIRS=(str(ppath('www')),)


MEDIA_BUNDLES= (
    ('main.css',
    'css/basic.css',
    'css/custom-theme/jquery-ui-1.9.2.custom.css',
    'css/style.css',
    'qtip2/jquery.qtip.css',
    ),
     
    ('print.css',
     'css/print.css'),
 
    ('jquery.js',
     'js/jquery-1.8.3.js',
     'js/jquery-ui-1.9.2.custom.min.js',
     'qtip2/jquery.qtip.js'
     ),
     ('processing.js',
      'js/processing.js',
     ),
    )
GENERATED_MEDIA_DIR = ppath('_generated_media',create=True)
GENERATED_MEDIA_NAMES_FILE = path(__file__).dirname() / '_generated_media_names.py'



# Make this unique, and don't share it with anybody.
SECRET_KEY = ')=m+&gn97_#2s!gi=ivgp%9j92cmj76+ecy&*kin#p&f%2p$78'

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
#    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
#     'django.template.loaders.eggs.Loader',
)

MIDDLEWARE_CLASSES = (
     # Media middleware has to come first
    'mediagenerator.middleware.MediaMiddleware',
   
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.auth.middleware.RemoteUserMiddleware',    
    'django.contrib.messages.middleware.MessageMiddleware',
    'npactweb.middleware.RedirectExceptionHandler',
)


TEMPLATE_CONTEXT_PROCESSORS = (
    #The defaults
    #"django.contrib.auth.context_processors.auth",
    "django.core.context_processors.debug",
    "django.core.context_processors.i18n",
    "django.core.context_processors.media",
    "django.core.context_processors.static",
    "django.contrib.messages.context_processors.messages",
    #addons
    "npactweb.context_processors.resolvermatch",

     )

ROOT_URLCONF = 'npactweb.urls'

#since we use the django.template.loaders.app_directories.Loader it automatically looks for apps' templates directory.
TEMPLATE_DIRS = ()

SESSION_ENGINE="django.contrib.sessions.backends.file"

INSTALLED_APPS = (
    #'django.contrib.auth',
    #'django.contrib.contenttypes',
    'django.contrib.sessions',
    #'django.contrib.sites',
    'django.contrib.messages',
    
    'mediagenerator',
    'npactweb'
    #'django.contrib.staticfiles',
    # Uncomment the next line to enable the admin:
    # 'django.contrib.admin',
    # Uncomment the next line to enable admin documentation:
    # 'django.contrib.admindocs',
)

DATABASES = {}

MESSAGE_STORAGE='django.contrib.messages.storage.cookie.CookieStorage'

#this will be passed to logging.dictConfig (or something equivalent)
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'verbose': {
            'format': '%(asctime)s %(process)d %(thread)d %(name)s %(levelname)s %(message)s'
            },
        'tty': {
            'format': "%(asctime)s.%(msecs)d %(name)-10s %(levelname)-8s %(message)s",
            'datefmt': '%H:%M:%S'
            },
        'simple': {
            'format': '%(levelname)s %(message)s'
            },
        'named_processes': {
            'format': "%(asctime)s %(process)s/%(name)s %(levelname)-8s %(message)s"
            },
    },
    'handlers': {
        'mail_admins': {
            'level': 'ERROR',
            'class': 'django.utils.log.AdminEmailHandler'
        },
        'console': {
            'level':'DEBUG',
            'class':'logging.StreamHandler',
            'formatter': 'tty',
        },
        'mainlogfile': {
            'class':'logging.handlers.WatchedFileHandler',
            'filename': ppath("logs",True) / "main.log",
            'formatter': 'verbose'
        },
        'django-tqdaemon': {
            'class':'logging.handlers.WatchedFileHandler',
            'filename': ppath("logs",True) / "django-tqdaemon.log",
            'formatter': 'named_processes',
        },
    },
    'loggers': {
        '': {
            'level': 'DEBUG',
            'handlers': ['mainlogfile'],
            },
        'django.request': {
            'handlers': ['mail_admins'],
            'level': 'ERROR',
            'propagate': True,
            },
        'cleanup': {
            'handlers': ['console'],
            'level': 'WARNING',
            'propagate': False,
        },
        # 'npact': {
        #     'handlers': ['console'],
        #     'level': 'DEBUG',
        #     'propagate': True,
        #     },
        # 'pynpact': {
        #     'handlers': ['console'],
        #     'level': 'DEBUG',
        #     'propagate': True,
        #     }
        'taskqueue': {
            'propagate': False,
            'level': 'DEBUG',
            'handlers': ['django-tqdaemon']
        },
        'multiprocessing': {
            'propagate': False,
            'level': 'INFO',
            'handlers': ['django-tqdaemon']
        }
    }
}

if DEBUG:
    def add_console_to(logger_name):
        LOGGING['loggers'][logger_name]['handlers'].append('console')
    add_console_to('')
    add_console_to('taskqueue')
    add_console_to('multiprocessing')
