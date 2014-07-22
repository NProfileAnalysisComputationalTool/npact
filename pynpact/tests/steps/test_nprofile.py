import os.path
import pytest
import py
from pynpact.steps import nprofile


def test_binfile_exists():
    assert nprofile.BIN
    assert os.path.exists(nprofile.BIN)


def test_plan(gbkconfig, plan_processor):
    "Test that the enqueue runs and returns"
    config = plan_processor(nprofile.plan, gbkconfig)
    filename = config[nprofile.OUTPUTKEY]
    assert filename
    p = py.path.local(filename)
    assert p.exists()
    # based on the default step:51, period:3 and the genome length of testgbk
    assert 59 == len(p.readlines())
