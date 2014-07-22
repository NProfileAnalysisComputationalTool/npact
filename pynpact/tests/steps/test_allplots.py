import os.path
import pytest
import py
import StringIO
from pynpact import parsing
from pynpact.steps import allplots


def test_binfile_exists():
    assert allplots.BIN
    assert os.path.exists(allplots.BIN)


def test_write_allplotsdef(gbkconfig):
    parsing.length(gbkconfig)
    buf = StringIO.StringIO()
    allplots.write_allplots_def(buf, gbkconfig, 1)
    lines = buf.getvalue().split('\n')
    assert 26 == len(lines)


def test_plan_allplots(gbkconfig, executor):
    allplots.allplots(gbkconfig, executor)
    assert len(gbkconfig['psnames']) > 0
    assert gbkconfig['psnames'][0] is not None


def test_combine_ps(gbkconfig, executor):
    allplots.combine_ps_files(gbkconfig, executor)
    assert gbkconfig['combined_ps_name']

