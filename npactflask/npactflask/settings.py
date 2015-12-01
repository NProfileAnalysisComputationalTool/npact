from path import Path
import os

DEBUG = os.environ.get('ENV', 'dev').lower() not in ('prod', 'production')

SEND_FILE_MAX_AGE_DEFAULT = 0 if DEBUG else 86400

MAIL_DEFAULT_SENDER = 'npact1.0@gmail.com'

WEBROOT = (Path(__file__).dirname() / "../../webroot").realpath()
LOGDIR = WEBROOT / 'logs'


def ppath(rel, create=True):
    abspath = (WEBROOT / rel).realpath()
    if abspath.exists():
        return abspath
    elif create:
        abspath.makedirs_p()
        return abspath
    else:
        raise Exception("Path '%s' doesn't exist." % abspath)


def get_secret():
    "Get the secret key for session signing"
    keyfile = WEBROOT / 'secret-key'
    if keyfile.exists():
        return keyfile.bytes()
    b = os.urandom(24)
    keyfile.write_bytes(b)
    return b

SECRET_KEY = get_secret()

UPLOADS = ppath('uploads')

# How many days should we keep upload files and products that haven't
# been accessed
ATIME_DEFAULT = 30
