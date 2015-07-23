import logging
import subprocess
from npactflask import app

from path import path

from pynpact import capproc
from npactflask.views import library_root

logger = app.logger

CLEANUP_PATHS = [path(app.config['UPLOADS']).realpath(),
                 path(app.config['TASKQUEUE']).realpath()]


def cleanup_old_files(days):
    days = int(days)
    err = False
    for p in CLEANUP_PATHS:
        err = err or 0 != clean_path(p, days)
    return not err


def clean_path(path, days):
    "This function actually runs the find and delete."
    logger.info("Cleaning older than: %d @ %r", days, path)
    cmd = ["find", path, "-atime", "+" + str(days), "-delete"]
    return capproc.capturedCall(cmd, logger=logger, stdin=False,
                                stderr_level=logging.WARNING)


def report_file_size():
    proc = capproc.capturedPopen(['du', '-h', '-s'] + CLEANUP_PATHS,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
    return proc.communicate()


def clear_library():
    logger.debug("Deleting %d from %s",
                 len(path(library_root()).files()), library_root())
    for f in path(library_root()).files():
        f.unlink()
    for d in path(library_root()).dirs():
        d.rmtree()
