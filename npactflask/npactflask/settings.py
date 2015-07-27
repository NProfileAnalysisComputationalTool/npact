from path import path

from npactflask import app

# TODO: Figure out mail configuration
EMAIL_USE_TLS = True
EMAIL_HOST = 'smtp.gmail.com'
EMAIL_PORT = 587
EMAIL_HOST_USER = 'npact1.0@gmail.com'
EMAIL_HOST_PASSWORD = 'sictransit2'


WEBROOT = (path(__file__).dirname() / "../../webroot").realpath()
if not WEBROOT.exists():
    raise Exception("Couldn't find webroot at %s" % WEBROOT)
else:
    app.config['WEBROOT'] = WEBROOT


def ppath(rel, create=True):
    abspath = (WEBROOT / rel).realpath()
    if abspath.exists():
        return abspath
    elif create:
        abspath.makedirs_p()
        return abspath
    else:
        raise Exception("Path '%s' doesn't exist." % abspath)

app.config['UPLOADS'] = ppath('uploads')


# How many days should we keep upload files and products that haven't
# been accessed
app.config['ATIME_DEFAULT'] = 30
