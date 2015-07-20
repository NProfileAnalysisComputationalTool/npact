DEBUG = True
import os
from path import path

WEBROOT = (path(__file__).dirname() / "../../webroot").realpath()
if not WEBROOT.exists():
    raise Exception("Couldn't find webroot at %s" % WEBROOT)


def ppath(rel, create=False):
    abspath = (WEBROOT / rel).realpath()
    if abspath.exists():
        return abspath
    elif create:
        os.makedirs(abspath)
        return abspath
    else:
        raise Exception("Path '%s' doesn't exist." % abspath)



# https://docs.djangoproject.com/en/dev/topics/email/#smtp-backend
EMAIL_USE_TLS = True
EMAIL_HOST = 'smtp.gmail.com'
EMAIL_PORT = 587
EMAIL_HOST_USER = 'npact1.0@gmail.com'
EMAIL_HOST_PASSWORD = 'sictransit2'


TIME_ZONE = 'America/New_York'


# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/home/media/media.lawrence.com/media/"
MEDIA_ROOT = ppath('uploads', True)

TQ_DIR = ppath('taskqueue', True)
import taskqueue
taskqueue.BASE_DIR = TQ_DIR


# how many days should we keep uploaded files and products that
# haven't been accessed before we delete them.
MEDIA_RETAIN_FOR = 7
ATIME_DEFAULT = 30


MIDDLEWARE_CLASSES = (
    # Media middleware has to come first
    'npactweb.middleware.NPactResponseExceptionHandler',
)


TEMPLATE_CONTEXT_PROCESSORS = (
    "npactweb.context_processors.resolvermatch",
    "npactweb.context_processors.addgetrelpath"
)


# this will be passed to logging.dictConfig (or something equivalent)
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
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'tty',
        },
        'mainlogfile': {
            'class': 'logging.handlers.WatchedFileHandler',
            'filename': ppath("logs", True) / "main.log",
            'formatter': 'verbose'
        },
        'django-tqdaemon': {
            'class': 'logging.handlers.WatchedFileHandler',
            'filename': ppath("logs", True) / "django-tqdaemon.log",
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
