import os.path
from pynpact import parsing
from pynpact.steps import allplots
from pynpact.util import which


def test_binfile_exists():
    assert allplots.BIN
    assert os.path.exists(allplots.BIN)


def test_build_allplotsdef(gbkconfig):
    parsing.length(gbkconfig)
    lines = allplots.build_allplots_def(gbkconfig, 1).split('\n')
    assert 26 == len(lines)


def test_plan_allplots(gbkconfig, executor):
    allplots.allplots(gbkconfig, executor)
    assert len(gbkconfig['psnames']) > 0
    assert gbkconfig['psnames'][0] is not None


def test_combine_ps(gbkconfig, executor):
    tasks = allplots.combine_ps_files(gbkconfig, executor)
    assert len(tasks) == 1
    psname = tasks[0]
    assert gbkconfig['combined_ps_name']
    assert gbkconfig['combined_ps_name'] == psname
    assert gbkconfig['allplots_result'] == gbkconfig['combined_ps_name']


def test_ps_to_pdf(gbkconfig, null_executor):
    jobs = allplots.convert_ps_to_pdf(gbkconfig, null_executor)
    assert len(jobs) == 1
    pdf_filename = jobs[0]
    assert pdf_filename == gbkconfig['pdf_filename']
    assert gbkconfig['pdf_filename'] == gbkconfig['allplots_result']


def test_all_the_way(gbkconfig, executor):
    filename = allplots.plan(gbkconfig, executor)[0]
    if which('ps2pdf'):
        assert filename == gbkconfig['pdf_filename']
    else:
        assert filename == gbkconfig['combined_ps_name']


def test_all_the_way_async(gbkconfig, async_executor):
    filename = allplots.plan(gbkconfig, async_executor)[0]
    if which('ps2pdf'):
        assert filename == gbkconfig['pdf_filename']
    else:
        assert filename == gbkconfig['combined_ps_name']
    async_executor.result(filename, timeout=10)
