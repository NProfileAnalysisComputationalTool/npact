PROC_TITLE = 'npact-taskqueue-daemon'

LISTEN_ADDRESS = ('127.0.0.1', 57129)
AUTH_KEY = 'npact'
BASE_DIR = '/tmp/'


import logging

from multiprocessing.managers import SyncManager
from functools import wraps


log = logging.getLogger(__name__)


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
    handler.setFormatter(logging.Formatter(
        "%(asctime)s %(processName)s/%(name)s %(levelname)-8s %(message)s",
        datefmt='%H:%M:%S'))
    tq_logger.addHandler(handler)
    mp_logger = multiprocessing.get_logger()
    mp_logger.addHandler(handler)
    mp_logger.setLevel(logging.INFO)


class ServerManager(SyncManager):
    pass

ServerManager.register('Server')


def memoize_noarg(fctn):
    "A custom memoizer for a no arg function; only need to store one value"
    @wraps(fctn)
    def memo():
        if memo.__cache__ is memo:
            memo.__cache__ = fctn()
        return memo.__cache__
    memo.__cache__ = memo
    return memo


@memoize_noarg
def instantiateit():
    from taskqueue.server import Server
    return Server()

# TODO: better name for logger param
def get_ServerManager(address=None, make_server=False, logger=False):
    if address is None:
        address = LISTEN_ADDRESS
    sm = ServerManager(address=address, authkey=AUTH_KEY)
    if logger:
        setup_logger(True)
    if make_server:
        method_to_typeid = {
            'get_task': 'AsyncResult'
        }
        try:
            # The callable in register is called whenever a client asks
            # for that, so need to instantiate it once and always return
            # that in order to have it actually shared.
            sm.register(
                'Server',
                callable=instantiateit, method_to_typeid=method_to_typeid)
            log.info("Server configured for socket at %s", sm.address)
        except:
            log.exception("Error registering taskqueue server")
    else:
        log.info("Returning a ServerManager client")
    return sm
