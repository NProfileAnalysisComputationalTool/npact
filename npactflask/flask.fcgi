#!/usr/bin/env python
#http://flask.pocoo.org/docs/0.10/deploying/fastcgi/#creating-a-fcgi-file
from flup.server.fcgi import WSGIServer
from npactflask import app

if __name__ == '__main__':
#TODO: Test and uncomment the following log code.
# from logging.handlers import WatchedFileHandler
# from logging import getLogger, Formatter

# vFormatter = Formatter('%(asctime)s %(process)d:%(thread)d %(name)-15s'
#                        ' %(levelname)-7s| %(message)s')

# root = getLogger('')
# roothandler = WatchedFileHandler(ppath('logs') / 'main.log')
# roothandler.setFormatter(vFormatter)
# root.addHandler(roothandler)

    #TODO: this probably needs the DispatcherMiddleware similar to runserver
    WSGIServer(app).run()
