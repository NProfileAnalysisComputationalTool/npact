#!/bin/sh

HOST="genome"
#make sure root doesn't have spaces in it or we'll have problems.
ROOT="/var/www/html/genome.ufl.edu/npact"

#HOST="progden"

rsync -vrp --exclude-from=.gitignore . $HOST:"${ROOT}"
ssh $HOST "[ -e ${ROOT}/webroot/uploads ] ||  mkdir -m 777 ${ROOT}/webroot/uploads"
ssh $HOST "[ -e ${ROOT}/webroot/logs ] || mkdir -m 777 ${ROOT}/webroot/logs"
ssh $HOST "${ROOT}/bootstrap.py"
ssh $HOST "env DJANGO_SETTINGS_MODULE='settings' ${ROOT}/ve/bin/django-admin.py generatemedia"
ssh $HOST "chmod -R 666 ${ROOT}/webroot/logs/*"
