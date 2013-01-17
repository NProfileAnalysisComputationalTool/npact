import os
import os.path
import logging
import signal
import subprocess
import re
import time

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
        except OSError:
            plf.break_lock()
            return False
    return pid

def stop():
    "Kill a running daemon instance"

    pid = status()
    if not pid:
        logger.warn("Daemon not running.")
        return INSTANCE_NOT_RUNNING

    # Try killing the daemon process
    try:
        for i in range(0, 30):
            os.kill(pid, signal.SIGTERM)
            time.sleep(0.5)
    except OSError, err:
        err = str(err)
        if err.find("No such process") > 0:
            get_pidfile().break_lock()
        else:
            logger.exception("Failed to stop daemon.")
            return OPERATION_FAILED
    return OPERATION_SUCCESSFUL

def restart():
    "Restart a daemon process"
    if status():
        kill_status = stop()
        if kill_status == OPERATION_FAILED:
            return kill_status
    return daemonize()

def kill(sig=signal.SIGKILL):
    uid = os.getuid()
    proc = subprocess.Popen(['ps', 'x', '-U', str(uid), '-o', 'pid,command'], stdout=subprocess.PIPE)
    lines = proc.stdout.readlines()[1:]
    logging.debug("Found %d processes for this user.", len(lines))
    for l in lines:
        l = l.strip()
        m = re.match('(\\d+) (npact-.*)', l)
        if m:
            pid,name = m.groups()
            logging.warning("Killing proc %s %r", pid, name)
            os.kill(int(pid), sig)



def daemonize():
    "Start a daemonized taskqueue"
    import daemon

    ##Before daemonizing figure out which handlers we want to keep.
    ##We only want to keep handlers for this library that aren't going
    ##to stderr or stdout
    other_handlers = set(logging._handlerList)
    fds = []
    l = logger
    while l:
        for h in l.handlers:
            if hasattr(h, 'stream') and \
              hasattr(h.stream, 'fileno') and \
              h.stream.fileno() not in [1,2]:
                fds.append(h.stream.fileno())
                other_handlers.discard(h)
        l = l.propagate and l.parent

    ## Kill of any other loggers
    logging.raiseExceptions = False
    logger.debug("Killing other loggers")
    logging.shutdown(list(other_handlers))
    #logger.debug("Remaining loggers: %r:%r", fds, logging._handlerList)

    pidfile = get_pidfile()
    if pidfile.is_locked():
        logger.error("Daemon already running")
        return INSTANCE_ALREADY_RUNNING
    else:
        logger.debug("Daemonizing, pidfile: %r", pidfile.path)
    try:
        with daemon.DaemonContext(pidfile=pidfile, files_preserve=fds, detach_process=True):
            logging.raiseExceptions = True
            try:
                import setproctitle
                import taskqueue
                setproctitle.setproctitle(taskqueue.PROC_TITLE)
                logger.debug("Successfully setproctitle: %r", taskqueue.PROC_TITLE)
            except:
                logger.exception("Couldn't setproctitle.")
                pass

            import taskqueue.server
            logger.info("Daemonized context")
            taskqueue.server.start_everything()
    finally:
        logger.warning("Exiting (hopefully intentionally)")


def start():
    "alias for daemonize"
    return daemonize()

def run():
    "Start the server in the foreground instead of as a daemon"
    import taskqueue.server
    logger.info("Running in foreground forever")
    taskqueue.server.start_everything()
