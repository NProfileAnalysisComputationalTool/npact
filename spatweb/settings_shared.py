# Django settings for spatweb project.
import os

DEBUG = True
TEMPLATE_DEBUG = DEBUG

from path import path

#physical path we're running this at, the idea is we can update only
#this variable on deploy and everything else should work out.
#PPATH="/home/ACCELERATION/nathan/projects/spat/spatweb/"
PPATH=path(__file__).dirname()
def ppath(rel, create=False) :
    abspath = (PPATH / rel).realpath()
    if abspath.exists():
        return abspath
    elif create:
        os.makedirs(abspath)
        return abspath
    else:
        raise Exception("Path '%s' doesn't exist." % abaspath)

from settings_logging import *

ADMINS = (
    ('Nathan Bird', 'nathan@acceleration.net'),
)

MANAGERS = ADMINS

DATABASES = {
    # 'default': {
    #     'ENGINE': 'django.db.backends.', # Add 'postgresql_psycopg2', 'postgresql', 'mysql', 'sqlite3' or 'oracle'.
    #     'NAME': 'spat',                      # Or path to database file if using sqlite3.
    #     'USER': 'spat',                      # Not used with sqlite3.
    #     'PASSWORD': 'Sp4t4g0r14l',                  # Not used with sqlite3.
    #     'HOST': '',                      # Set to empty string for localhost. Not used with sqlite3.
    #     'PORT': '',                      # Set to empty string for default. Not used with sqlite3.
    # }
}

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# On Unix systems, a value of None will cause Django to use the same
# timezone as the operating system.
# If running in a Windows environment this must be set to the same as your
# system time zone.
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

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash.
# Examples: "http://media.lawrence.com/media/", "http://example.com/media/"
MEDIA_URL = ''


#### django-mediagenerator settings

MEDIA_DEV_MODE = DEBUG
DEV_MEDIA_URL = '/devassets/'
PRODUCTION_MEDIA_URL = '/assets/'
GLOBAL_MEDIA_DIRS=(str(ppath('www')),)


MEDIA_BUNDLES= (
    ('main.css',
     'css/basic.css',
     'css/style.css',
     'css/custom-theme/jquery-ui-1.8.16.custom.css',
     'qtip2/jquery.qtip.css',
    ),

    ('print.css',
     'css/print.css'),

    ('jquery.js',
     'js/jquery-1.6.4.min.js',
     'js/jquery-ui-1.8.16.custom.min.js',
     'qtip2/jquery.qtip.js',),
    )




# Make this unique, and don't share it with anybody.
SECRET_KEY = ')=m+&gn97_#2s!gi=ivgp%9j92cmj76+ecy&*kin#p&f%2p$78'

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
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
    'django.contrib.messages.middleware.MessageMiddleware',
    'spat.middleware.RedirectExceptionHandler',
)

#What it is by default.
# TEMPLATE_CONTEXT_PROCESSORS = (
#     "django.contrib.auth.context_processors.auth",
#     "django.core.context_processors.debug",
#     "django.core.context_processors.i18n",
#     "django.core.context_processors.media",
#     "django.core.context_processors.static",
#     "django.contrib.messages.context_processors.messages")
ROOT_URLCONF = 'spatweb.urls'

TEMPLATE_DIRS = (
    ppath('templates'),

    # Put strings here, like "/home/html/django_templates" or "C:/www/django/templates".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
)

SESSION_ENGINE="django.contrib.sessions.backends.file"

INSTALLED_APPS = (
    #'django.contrib.auth',
    #'django.contrib.contenttypes',
    'django.contrib.sessions',
    #'django.contrib.sites',
    'django.contrib.messages',
    
    'mediagenerator',
    'spat'
    #'django.contrib.staticfiles',
    # Uncomment the next line to enable the admin:
    # 'django.contrib.admin',
    # Uncomment the next line to enable admin documentation:
    # 'django.contrib.admindocs',
)


MESSAGE_STORAGE='django.contrib.messages.storage.cookie.CookieStorage'

EXC_TRACE_PATH=ppath("logs/exceptions", True)
