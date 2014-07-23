import os.path
import pytest
import py
from pynpact.steps import acgt_gamma


def test_binfile_exists():
    assert acgt_gamma.BIN
    assert os.path.exists(acgt_gamma.BIN)


def test__acgt_gamma(gbkconfig, tmpdir):
    outputdir = tmpdir.join('predict')
    acgt_gamma._acgt_gamma(gbkconfig, str(outputdir))
    assert 0 < len(outputdir.listdir())


def test_plan(gbkconfig, executor):
    acgt_gamma.plan(gbkconfig, executor)

    predictdir = gbkconfig[acgt_gamma.OUTPUTKEY]
    assert predictdir
    p = py.path.local(predictdir)
    assert p.exists()


def test_plan_async(gbkconfig, async_executor):
    acgt_gamma.plan(gbkconfig, async_executor)

    predictdir = gbkconfig[acgt_gamma.OUTPUTKEY]

    assert predictdir
    assert predictdir == async_executor.result(predictdir, 1)
    p = py.path.local(predictdir)
    assert p.exists()
