import os.path
import pytest
import py
import StringIO
from pynpact.steps import allplots


def test_binfile_exists():
    assert allplots.BIN
    assert os.path.exists(allplots.BIN)


@pytest.fixture()
def apstep(gbkconfig, tmpdir, executor):
    return allplots.AllplotsStep(
        config=gbkconfig, outputdir=str(tmpdir), executor=executor)


def test_apstep(apstep):
    assert apstep

def test_write_allplotsdef(gbkconfig):
    buf = StringIO.StringIO()
    allplots.write_allplots_def(buf, gbkconfig, 1)
    lines =  buf.getvalue().split('\n')
    assert 26 == len(lines)


def test_enqueue(apstep):
    filenames = apstep.enqueue()
    assert filenames
    assert 1 == len(filenames)
    p = py.path.local(filenames[0])
    assert p.exists()
