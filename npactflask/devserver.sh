#!/usr/bin/env sh

export ENV=dev
gunicorn --bind=127.0.0.1:5000 --reload --worker-class=gevent --access-logfile=- server
