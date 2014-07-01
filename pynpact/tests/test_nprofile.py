import os.path
import pytest
import py
from pynpact.steps import nprofile


def test_binfile_exists():
    assert nprofile.BIN
    assert os.path.exists(nprofile.BIN)


@pytest.fixture()
def npstep(gbkconfig, tmpdir, executor):
    return nprofile.NprofileStep(
        config=gbkconfig, outputdir=str(tmpdir), executor=executor)


def test_npstep(npstep):
    assert npstep


def test_enqueue(npstep):
    "Test that the enqueue runs and returns"
    filename = npstep.enqueue()
    assert filename
    p = py.path.local(filename)
    assert p.exists()
    #based on the default step:51, period:3 and the genome length of testgbk
    assert 59 == len(p.readlines())
