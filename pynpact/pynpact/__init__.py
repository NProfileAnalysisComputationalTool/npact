import logging
import os.path

from pkg_resources import resource_filename

class NullHandler(logging.Handler):
    def emit(self, record):
        pass

logging.getLogger('pynpact').addHandler(NullHandler())


BINPATH=os.path.realpath(os.path.join(os.path.dirname(__file__),"bin"))

def binfile(name) :
    path = os.path.join(BINPATH, name)
    assert os.path.exists(path), "Missing pynpact/bin/" + name + ", perhaps you need to run make?"
    return  path


DATAPATH=resource_filename(__name__, 'data')
