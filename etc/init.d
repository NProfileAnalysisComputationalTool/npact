#!/bin/sh
#
# /etc/rc.d/init.d/npact-supervisor
#
# Supervisor is a client/server system that
# allows its users to monitor and control a
# number of processes on UNIX-like operating
# systems.
#
# chkconfig: - 64 36
# description: NPACT Supervisor process controller
# processname: supervisord
# pidfile: /var/run/npact-supervisord.pid

### BEGIN INIT INFO
# Provides: npact-supervisor
# Required-Start: $local_fs $network
# Required-Stop: $local_fs $network
# Default-Start:  2 3 4 5
# Default-Stop: 0 1 6
# Short-Description: NPACT Supervisor process controller
# Description: NPACT Supervisor process controller
### END INIT INFO


# Source init functions
. /etc/rc.d/init.d/functions


RETVAL=0
prog="NPACT supervisord"
PREFIX="/var/www/html/genome.ufl.edu/npact"
SUPERVISORD="${PREFIX}/ve/bin/supervisord"
PID_FILE="/var/run/npact-supervisord.pid"
CONFIG_FILE="${PREFIX}/etc/supervisord.conf"

start()
{
        echo -n $"Starting $prog: "
        $SUPERVISORD -c $CONFIG_FILE --pidfile $PID_FILE && success || failure
        RETVAL=$?
        echo
        return $RETVAL
}

stop()
{
        echo -n $"Stopping $prog: "
        killproc -p $PID_FILE -d 10 $SUPERVISORD
        RETVAL=$?
        echo
}

reload()
{
        echo -n $"Reloading $prog: "
        if [ -n "`pidfileofproc $SUPERVISORD`" ] ; then
            killproc $SUPERVISORD -HUP
        else
            # Fails if the pid file does not exist BEFORE the reload
            failure $"Reloading $prog"
        fi
        sleep 1
        if [ ! -e $PID_FILE ] ; then
            # Fails if the pid file does not exist AFTER the reload
            failure $"Reloading $prog"
        fi
        RETVAL=$?
        echo
}

case "$1" in
        start)
                start
                ;;
        stop)
                stop
                ;;
        restart)
                stop
                start
                ;;
        reload)
                reload
                ;;
        status)
                status -p $PID_FILE $SUPERVISORD
                RETVAL=$?
                ;;
        *)
                echo $"Usage: $0 {start|stop|restart|reload|status}"
                RETVAL=1
esac
exit $RETVAL
