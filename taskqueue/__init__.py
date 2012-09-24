import logging

PERSISTENCE_DIR = '/tmp/taskqueue'

LISTEN_ADDRESS = ('localhost', 57129)

AUTH_KEY = 'npact'


class NoSuchTaskError(Exception):
    pass


def setup_logger(verbose):
    tq_logger = logging.getLogger('taskqueue')
    tq_logger.propagate=False
    tq_logger.setLevel(logging.DEBUG if verbose else logging.WARNING)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter("%(asctime)s %(name)-10s "
                                           "%(levelname)-8s %(message)s",
                                           datefmt='%H:%M:%S'))
    tq_logger.addHandler(handler)
