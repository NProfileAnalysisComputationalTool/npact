#!/bin/sh

HOST="lucianob@genome"
#make sure root doesn't have spaces in it or we'll have problems.
ROOT="/var/www/html/genome.ufl.edu/npact"


rsync -vrt --exclude-from=.gitignore . $HOST:"${ROOT}"
ssh $HOST "[ -e ${ROOT}/webroot/uploads ] ||  mkdir -p -m 777 ${ROOT}/webroot/uploads"
ssh $HOST "[ -e ${ROOT}/webroot/logs ] || mkdir -p -m 777 ${ROOT}/webroot/logs"
ssh $HOST "umask 0002; ${ROOT}/bootstrap.py"
ssh $HOST "umask 0002; env DJANGO_SETTINGS_MODULE='settings' ${ROOT}/ve/bin/django-admin.py generatemedia"
ssh $HOST "umask 0002; chgrp -R webadmin $ROOT; chmod g+w -R $ROOT"
ssh $HOST "chmod -R a+w ${ROOT}/webroot"
echo "===== INFO ====="
ssh $HOST "[ -e ${ROOT}/.htpasswd ]" || echo "WARN: Missing .htpasswd file; please create."
echo "INFO: Go to http://genome.ufl.edu/npact/management and restart the daemon."
