"Provide a bunch of test fixtures for pynpact testing"
import logging
import pytest
from path import path

from pynpact import parsing, executors


@pytest.fixture(scope="module", autouse=True)
def setup_logging():
    logging.root.setLevel(logging.DEBUG)
    logging.basicConfig()


TESTGBK = "testdata/NC_017123.gbk"
TESTFNA = "testdata/NC_017123.fna"


@pytest.fixture()
def executor():
    return executors.InlineExecutor()


@pytest.fixture()
def gbkfile():
    gbk = path(__file__).dirname().joinpath(TESTGBK)
    assert gbk.exists()
    return str(gbk)


@pytest.fixture()
def gbkconfig(gbkfile, tmpdir):
    return parsing.initial(gbkfile, outputdir=str(tmpdir))


@pytest.fixture()
def fnafile():
    fna = path(__file__).dirname().joinpath(TESTFNA)
    assert fna.exists()
    return str(fna)


@pytest.fixture()
def fnaconfig(fnafile, tmpdir):
    return parsing.initial(fnafile, outputdir=str(tmpdir))
