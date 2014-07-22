import os.path
import pytest
import py
from pynpact.steps import extract


def test_binfile_exists():
    assert extract.BIN
    assert os.path.exists(extract.BIN)


def test_plan(gbkconfig, plan_processor):
    config = plan_processor(extract.plan, gbkconfig)

    filename = config['File_of_published_accepted_CDSs']
    assert filename
    p = py.path.local(filename)
    assert p.exists()
    # based on how many genes are in testgbk
    assert 3 == len(p.readlines())
