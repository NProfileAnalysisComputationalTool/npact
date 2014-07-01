import os.path
import pytest
import py
from pynpact.steps import extract


def test_binfile_exists():
    assert extract.BIN
    assert os.path.exists(extract.BIN)


@pytest.fixture()
def exstep(gbkconfig, tmpdir, executor):
    return extract.ExtractStep(
        config=gbkconfig, outputdir=str(tmpdir), executor=executor)


def test_exstep(exstep):
    assert exstep


def test_enqueue(exstep):
    filename = exstep.enqueue()
    assert filename
    p = py.path.local(filename)
    assert p.exists()
    # based on how many genes are in testgbk
    assert 3 == len(p.readlines())
