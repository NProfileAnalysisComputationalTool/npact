
# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'verbose': {
            'format': '%(asctime)s %(process)d %(thread)d %(module)s %(levelname)s %(message)s'
            },
        'tty': {
            'format': "%(asctime)s %(module)-10s %(levelname)-8s %(message)s",
            'datefmt': '%H:%M:%S'
            },
        'simple': {
            'format': '%(levelname)s %(message)s'
            },
    },
    'handlers': {
        'mail_admins': {
            'level': 'ERROR',
            'class': 'django.utils.log.AdminEmailHandler'
        },
        'console':{
            'level':'DEBUG',
            'class':'logging.StreamHandler',
            'formatter': 'tty',
        },
    },
    'loggers': {
        '': {
            'handlers':['console'],
            'level': 'DEBUG',
            },
        'django.request': {
            'handlers': ['mail_admins'],
            'level': 'ERROR',
            'propagate': True,
            },
        # 'spat': {
        #     'handlers': ['console'],
        #     'level': 'DEBUG',
        #     'propagate': True,
        #     },
        # 'pynpact': {
        #     'handlers': ['console'],
        #     'level': 'DEBUG',
        #     'propagate': True,
        #     }
        }
    }

EXC_TRACE_FNAME="%s"
