from path import path

DEBUG = False

# TODO: Figure out mail configuration
EMAIL_USE_TLS = True
EMAIL_HOST = 'smtp.gmail.com'
EMAIL_PORT = 587
EMAIL_HOST_USER = 'npact1.0@gmail.com'
EMAIL_HOST_PASSWORD = 'sictransit2'

SECRET_KEY = 'cf0cb53d-1ff1-4074-a65d-977831de66af'

WEBROOT = (path(__file__).dirname() / "../../webroot").realpath()


def ppath(rel, create=True):
    abspath = (WEBROOT / rel).realpath()
    if abspath.exists():
        return abspath
    elif create:
        abspath.makedirs_p()
        return abspath
    else:
        raise Exception("Path '%s' doesn't exist." % abspath)

UPLOADS = ppath('uploads')

# How many days should we keep upload files and products that haven't
# been accessed
ATIME_DEFAULT = 30
