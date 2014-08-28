import os.path
import pytest
import py
from pynpact.steps import extract


def test_binfile_exists():
    assert extract.BIN
    assert os.path.exists(extract.BIN)


def test_plan(gbkconfig, executor):
    extract.plan(gbkconfig, executor)

    filename = gbkconfig[extract.OUTPUTKEY]
    assert filename
    p = py.path.local(filename)
    assert p.exists()
    # based on how many genes are in testgbk
    assert 3 == len(p.readlines())


def test_plan_async(gbkconfig, async_executor):
    extract.plan(gbkconfig, async_executor)

    filename = gbkconfig[extract.OUTPUTKEY]
    assert filename
    async_executor.result(filename, 1)
    p = py.path.local(filename)
    assert p.exists()
    # based on how many genes are in testgbk
    assert 3 == len(p.readlines())
