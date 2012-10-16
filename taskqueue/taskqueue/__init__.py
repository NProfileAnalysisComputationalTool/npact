
PROC_TITLE = 'npact-taskqueue-daemon'

LISTEN_ADDRESS = ('localhost', 57129)
AUTH_KEY = 'npact'
BASE_DIR = '/tmp/'


class NoSuchTaskError(Exception):
    pass


def setup_logger(verbose):
    import logging
    import multiprocessing

    tq_logger = logging.getLogger('taskqueue')
    tq_logger.propagate = False
    tq_logger.setLevel(logging.DEBUG if verbose else logging.WARNING)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter("%(asctime)s %(processName)s/%(name)s "
                                           "%(levelname)-8s %(message)s",
                                           datefmt='%H:%M:%S'))
    tq_logger.addHandler(handler)
    mp_logger = multiprocessing.get_logger()
    mp_logger.addHandler(handler)
    mp_logger.setLevel(logging.INFO)

