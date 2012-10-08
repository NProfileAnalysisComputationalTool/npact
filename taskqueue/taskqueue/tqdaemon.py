import os
import os.path
import logging 

import lockfile
from lockfile import pidlockfile

import taskqueue

OPERATION_SUCCESSFUL = 0
OPERATION_FAILED = 1
INSTANCE_ALREADY_RUNNING = 2
INSTANCE_NOT_RUNNING = 3
SET_USER_FAILED = 4

logger = logging.getLogger('taskqueue.daemon')

def get_pidfile():
    return pidlockfile.PIDLockFile(
        os.path.join(taskqueue.BASE_DIR, 'tqdaemon.pid'))

def status():
    "Check whether daemon is already running. return pid or False"
    plf = get_pidfile()
    pid = plf.read_pid()
    if pid:
        try:
            os.kill(pid, 0) #no-op signal, only indicates whether it could be sent
        except OSError, e:
            plf.break_lock()
            return False
    return pid

def stop():
    "Kill a running daemon instance"

    pid = status()
    if not pid:
        return INSTANCE_NOT_RUNNING

    # Try killing the daemon process	
    try:
        for i in range(0, 30):
            os.kill(pid, SIGTERM)
            time.sleep(0.5)
    except OSError, err:
        err = str(err)
        if err.find("No such process") > 0:
            self.pidfile.break_lock()
        else:
            return OPERATION_FAILED
    return OPERATION_SUCCESSFUL

def restart():
    "Restart a daemon process"
    if status():
        kill_status = kill()
        if kill_status == OPERATION_FAILED:
            return kill_status
    return daemonize()


def daemonize():
    "Start a deaemonized taskqueue"
    import logging
    
    import daemon
    with daemon.DaemonContext(pidfile=get_pidfile()):
        import taskqueue.server
        logger.info("Daemonized context")
        taskqueue.server.start_everything()

def start():
    "alias for daemonize"
    return daemonize()

def run():
    "Start the server in the foreground instead of as a daemon"
    import taskqueue.server
    logger.info("Running in foreground forever")
    taskqueue.server.start_everything()
