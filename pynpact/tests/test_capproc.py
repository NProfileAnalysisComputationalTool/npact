import pytest
import cPickle

from pynpact import capproc


def test_calledProcessError():
    cpe = capproc.calledProcessError(11, 'cmd')
    s = cPickle.dumps(cpe)
    assert s
    cpeup = cPickle.loads(s)
    assert cpeup
    assert cpeup.returncode == cpe.returncode
    assert cpeup.cmd == cpe.cmd
