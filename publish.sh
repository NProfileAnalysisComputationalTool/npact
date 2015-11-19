#!/bin/sh

HOST="genome.ufl.edu"
#make sure root doesn't have spaces in it or we'll have problems.
ROOT="/var/www/html/genome.ufl.edu/npact"


ssh -A $HOST "cd $ROOT; umask 0002; git pull"
ssh $HOST "[ -e ${ROOT}/webroot/uploads ] ||  mkdir -p -m 777 ${ROOT}/webroot/uploads"
ssh $HOST "[ -e ${ROOT}/webroot/logs ] || mkdir -p -m 777 ${ROOT}/webroot/logs"
ssh $HOST "umask 0002; ${ROOT}/bootstrap.py"
ssh $HOST "chgrp --quiet -R webadmin $ROOT; chmod --quiet g+w -R $ROOT"
ssh $HOST "chmod --quiet a+w ${ROOT}/webroot"
echo "===== INFO ====="
ssh $HOST "[ -e ${ROOT}/.htpasswd ]" || echo "WARN: Missing .htpasswd file; please create."
cd $ROOT
ve/bin/supervisorctl status npactserver
