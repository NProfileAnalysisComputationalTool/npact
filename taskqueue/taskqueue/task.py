import cPickle
import logging
import errno
import sys
import os.path
import string
import random

from multiprocessing.managers import RemoteError
from multiprocessing import TimeoutError
from path import path

from taskqueue import NoSuchTaskError
from taskqueue.client import get_server

log = logging.getLogger(__name__)


def sanitize_id(tid):
    assert isinstance(tid, basestring)
    return tid.replace(os.path.sep, '_')


def randomid():
    length = 16
    return ''.join(
        random.choice(string.lowercase) for i in range(length))


class Task(object):
    func = None
    tid = None
    path = None

    def __repr__(self):
        return "#<Task(tid=%r)>" % (self.tid)

    def __init__(self, func, tid=None, after=None):
        self.func = func
        self.tid = tid
        self.after = after

    def run(self):
        if self.after:
            log.debug("Waiting on %d tasks.", len(self.after))
            server = get_server()
            tid = None
            try:
                while self.after:
                    tid = self.after.pop()
                    try:
                        server.get_task(tid).wait(1)
                    except TimeoutError:
                        self.after.append(tid)
                        log.debug("Task %r not ready")
            except RemoteError:
                log.info("Task %r erred, aborting", tid)
                raise

        log.debug("Task running.")
        return self.func()

    def pickle(self, work_dir):
        "Writes the task def to a file via pickle in case of server restart"

        work_dir = path(work_dir)
        while True:
            fd = None
            # Unless given a specific tid generate one randomly
            tid = self.tid or randomid()
            p = work_dir / sanitize_id(tid)
            log.debug("Pickling %r to %r", self, p)
            try:
                # write out the file, make sure we have unique handle to it
                fd = os.open(p, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                fd = os.fdopen(fd, 'w')
                cPickle.dump(self.func, fd, cPickle.HIGHEST_PROTOCOL)
                self.tid = tid
            except OSError as e:
                if e.errno == errno.EEXIST:
                    # file wasn't unique; if it was a specified id
                    # then we're already done, if was a random id try
                    # again
                    if tid is self.tid:
                        return self.tid
                    else:
                        log.debug("id collision; retrying")
                        continue
                else:
                    raise
            finally:
                if fd:
                    fd.close()

        return self.tid

    @classmethod
    def unpickle(klass, tid, work_dir):
        "Read a task definition from file."
        path = os.path.join(work_dir, sanitize_id(tid))
        if not os.path.exists(path):
            raise NoSuchTaskError()
        with open(path, 'r') as f:
            func = cPickle.load(f)
        return klass(func=func, tid=tid)


def async_wrapper(task, work_dir):
    """Wrapper func around executing a task; ensures env setup and teardown.

    * This is run inside a pool process via the pool apply_async above.
    * Ensures the job 'todo' file is marked finished
    * sets up logging to a file for this task
    """
    work_dir = path(work_dir)

    # Logging setup:
    # In this proc all logging should go to a file named after the task
    root_logger = logging.root
    # blow away any other handlers from previous uses of this process.
    root_logger.handlers = []
    log_path = work_dir / (task.tid + '.log')
    file_handler = logging.FileHandler(log_path, mode='a')
    file_handler.setFormatter(
        logging.Formatter(
            '%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
            datefmt='%Y%m%d %H:%M:%S'))

    root_logger.setLevel(logging.DEBUG)
    root_logger.addHandler(file_handler)

    stderrlog = work_dir / (sanitize_id(task.tid) + '.stderr')
    oldstderr = sys.stderr
    try:
        with open(stderrlog, 'w') as f:
            sys.stderr = f
            rc = task.run()
            log.info("Finished %r", task.tid)
            return rc
    finally:
        sys.stderr = oldstderr
