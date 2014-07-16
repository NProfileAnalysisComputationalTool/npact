"Provide a bunch of test fixtures for pynpact testing"
import random
import string

import pytest
from path import path

from pynpact import prepare

TESTGBK = "testdata/NC_017123.gbk"
TESTFNA = "testdata/NC_017123.fna"


def randomid(length=16):
    return ''.join(
        random.choice(string.lowercase) for i in range(length))


class InlineExecutor(object):
    "Does everything in line"
    pool = None
    tasks = None

    def __init__(self):
        self.tasks = {}

    def enqueue(self, callable, tid=None, after=None):
        if tid is None:
            tid = randomid()

        if after is not None:
            for aid in after:
                assert aid in self.tasks, \
                    "The InlineExecutor can't be after a task that doesn't exist yet"

        if tid not in self.tasks:
            self.tasks[tid] = callable()
        return tid

    def result(self, tid):
        return self.tasks[tid]


@pytest.fixture()
def executor():
    return InlineExecutor()


@pytest.fixture()
def gbkfile():
    gbk = path(__file__).dirname().joinpath(TESTGBK)
    assert gbk.exists()
    return str(gbk)


@pytest.fixture()
def gbkconfig(gbkfile):
    return prepare.default_config(gbkfile)


@pytest.fixture()
def fnafile():
    fna = path(__file__).dirname().joinpath(TESTFNA)
    assert fna.exists()
    return str(fna)


@pytest.fixture()
def fnaconfig(fnafile):
    return prepare.default_config(fnafile)
