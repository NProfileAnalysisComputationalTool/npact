import os
import os.path
import logging
import signal
import subprocess
import re
import time

from lockfile import pidlockfile

import taskqueue

OPERATION_SUCCESSFUL = 0
OPERATION_FAILED = 1
INSTANCE_ALREADY_RUNNING = 2
INSTANCE_NOT_RUNNING = 3
SET_USER_FAILED = 4

def get_pidfile():
    return pidlockfile.PIDLockFile(
        os.path.join(taskqueue.BASE_DIR, 'tqdaemon.pid'))


def status():
    "Check whether daemon is already running. return pid or False"
    plf = get_pidfile()
    pid = plf.read_pid()
    if pid and check_for_process(pid):
        return pid
    return False


def check_for_process(pid):
    if pid:
        try:
            # no-op signal, only indicates whether it could be sent
            os.kill(pid, 0)
            return True
        except OSError:
            get_pidfile().break_lock()
            return False


def stop():
    "Kill a running daemon instance"

    pid = status()
    if not pid:
        logging.getLogger('taskqueue.daemon').warn("Daemon not running.")
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
            logging.getLogger('taskqueue.daemon').exception("Failed to stop daemon.")
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
    proc = subprocess.Popen(
        ['ps', 'x', '-U', str(uid), '-o', 'pid,command'],
        stdout=subprocess.PIPE)
    lines = proc.stdout.readlines()[1:]
    logging.getLogger('taskqueue.daemon').debug(
        "Searching %d processes for this user.", len(lines))
    killed = 0
    for l in lines:
        l = l.strip()
        m = re.match('(\\d+) (npact-.*)', l)
        if m:
            pid, name = m.groups()
            logging.getLogger('taskqueue.daemon').warning(
                "Killing proc %s %r", pid, name)
            try:
                os.kill(int(pid), sig)
                killed += 1
            except:
                logging.getLogger('taskqueue.daemon').exception(
                    "Error killing %s %r", pid, name)
    return killed


def daemonize():
    """Start a daemonized taskqueue

    This the code run after the initial fork to help get us all the
    way to daemonized.

    """
    import daemon
    import taskqueue
    import logging
    logger = logging.getLogger('taskqueue.daemon')
    try:
        # make tqdaemon a bit nicer than whatever parent launched us.
        os.nice(4)
    except:
        pass
    # Want to completely stop logging in this proc
    pidfile = get_pidfile()
    check_for_process(pidfile.read_pid())
    if pidfile.is_locked():
        logger.error("Daemon already running")
        return INSTANCE_ALREADY_RUNNING
    else:
        logger.info("Daemonizing, pidfile: %r", pidfile.path)
    logging.shutdown()
    with daemon.DaemonContext(pidfile=pidfile, detach_process=True,
                              working_directory=taskqueue.BASE_DIR):
        tqdaemonlog()
        _proctitle()
        log = logging.getLogger('taskqueue.daemon')
        import taskqueue.server
        log.info("Daemonized context")
        sm = taskqueue.get_ServerManager(make_server=True)
        sm.get_server().serve_forever()


def _proctitle():
    import taskqueue
    import logging
    log = logging.getLogger('taskqueue.daemon')
    try:
        import setproctitle
        setproctitle.setproctitle(taskqueue.PROC_TITLE)
        log.debug("Successfully setproctitle: %r", taskqueue.PROC_TITLE)
    except:
        log.exception("Couldn't setproctitle.")
        pass


def tqdaemonlog():
    import logging
    import multiprocessing
    from logging.handlers import WatchedFileHandler
    logging.raiseExceptions = True
    handler = WatchedFileHandler(
        os.path.join(taskqueue.BASE_DIR, 'tqdaemon.log'))
    handler.setFormatter(logging.Formatter(
        "%(asctime)s %(processName)15s/%(module)-8s %(levelname)-8s %(message)s",
        datefmt='%H:%M:%S'))
    tqlog = logging.getLogger('taskqueue')
    mplog = multiprocessing.get_logger()
    tqlog.setLevel(logging.DEBUG)
    tqlog.propagate = False
    mplog.propagate = False
    tqlog.handlers = mplog.handlers = [handler]
    tqlog.info('Finished setting up logging')
    return tqlog


def start():
    "alias for daemonize"
    return daemonize()


def run():
    "Start the server in the foreground instead of as a daemon"
    import taskqueue.server
    logger.info("Running in foreground forever")
    taskqueue.get_ServerManager(make_server=True).get_server().serve_forever()
