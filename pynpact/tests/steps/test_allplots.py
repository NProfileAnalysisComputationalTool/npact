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


def test_plan_allplots(gbkconfig, plan_processor):
    config = plan_processor(allplots.allplots, gbkconfig)
    assert len(config['psnames']) > 0
    assert config['psnames'][0] is not None



def test_combine_ps(gbkconfig, plan_processor):
    config = plan_processor(allplots.combine_ps_files, gbkconfig)
    assert config['combined_ps_name']
