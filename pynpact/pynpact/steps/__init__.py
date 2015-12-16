import logging
from functools import wraps
from pynpact.util import mkstemp_rename, delay
from path import Path


log = logging.getLogger('pynpact.steps')


def enqueue(func, executor, config, target, after=None):
    """Enqueue the function on the executor if target doesn't exist.

    The function is expected to have been decorated by `@producer` and
    should have public interface of `(config, target_file)`.

    This function checks that the target file doesn't already exist
    and enqueues it with the correct parameters returning back the job
    id list.

    """
    return [
        executor.enqueue(
            delay(func)(config, target), tid=target, after=after)]


def producer(tmpmanager=mkstemp_rename):
    """A decorator for the functions that actually produce output

    These functions have public interface of (config, target_path) but
    are given parameters of (config, tmp_path). The
    tmp_path/target_path is managed by `mkstemp_rename`, or
    `mkdtemp_rename` (parameter to this decorator).
    """
    def getfn(func):
        @wraps(func)
        def wrapper(config, target):
            if not Path(target).exists():
                with tmpmanager(target, log=log) as tmp:
                    func(config, tmp)
            return target
        return wrapper
    return getfn
