import logging
import os.path

class NullHandler(logging.Handler):
    def emit(self, record):
        pass

logging.getLogger('pynpact').addHandler(NullHandler())


BINPATH=os.path.realpath(os.path.join(os.path.dirname(__file__),"bin"))

def binfile(name) :
    path = os.path.join(BINPATH, name)
    assert os.path.exists(path), "Missing pynpact/bin/" + name + ", perhaps you need to run make?"
    return  path


DATAPATH=os.path.realpath(os.path.join(os.path.dirname(__file__), "data"))

