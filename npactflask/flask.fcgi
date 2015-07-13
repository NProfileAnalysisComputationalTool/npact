#!/usr/bin/env python
#http://flask.pocoo.org/docs/0.10/deploying/fastcgi/#creating-a-fcgi-file
from flup.server.fcgi import WSGIServer
from npactflask import app

if __name__ == '__main__':
    WSGIServer(app).run()
