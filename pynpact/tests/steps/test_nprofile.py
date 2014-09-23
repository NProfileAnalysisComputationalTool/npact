import os.path
import pytest
import py
from pynpact.steps import nprofile


def test_binfile_exists():
    assert nprofile.BIN
    assert os.path.exists(nprofile.BIN)


def test_gbk(gbkconfig, executor):
    "Test that the enqueue runs and returns"
    nprofile.plan(gbkconfig, executor)
    filename = gbkconfig[nprofile.OUTPUTKEY]
    assert filename
    p = py.path.local(filename)
    assert p.exists()
    # based on the default step:51, period:3 and the genome length of testgbk
    assert 59 == len(p.readlines())


#@pytest.mark.xfail()
def test_fna(fnaconfig, executor):
    nprofile.plan(fnaconfig, executor)
    filename = fnaconfig[nprofile.OUTPUTKEY]
    assert filename
    p = py.path.local(filename)
    assert p.exists()
    # based on the default step:51, period:3 and the genome length of testgbk
    assert 59 == len(p.readlines())
