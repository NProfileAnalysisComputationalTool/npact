import os.path
import pytest
import py
import StringIO
from pynpact import parsing
from pynpact.steps import allplots
from pynpact.util import which


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
    psname = allplots.combine_ps_files(gbkconfig, executor)
    assert gbkconfig['combined_ps_name']
    assert gbkconfig['combined_ps_name'] == psname


def test__ps2pdf(patcher, tmpdir):
    mock = patcher.patch('pynpact.capproc.capturedCall')
    allplots._ps2pdf('foo', str(tmpdir.join('foo')))
    assert mock.called


def test_ps_to_pdf(gbkconfig, null_executor, patcher):
    pdf_filename = allplots.convert_ps_to_pdf(gbkconfig, null_executor)
    assert pdf_filename == gbkconfig['pdf_filename']


def test_all_the_way(gbkconfig, executor):
    filename = allplots.plan(gbkconfig, executor)
    if which('ps2pdf'):
        assert filename == gbkconfig['pdf_filename']
    else:
        assert filename == gbkconfig['combined_ps_name']


def test_all_the_way_async(gbkconfig, async_executor):
    filename = allplots.plan(gbkconfig, async_executor)
    if which('ps2pdf'):
        assert filename == gbkconfig['pdf_filename']
    else:
        assert filename == gbkconfig['combined_ps_name']
    async_executor.result(filename, timeout=10)
