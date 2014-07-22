import pytest

from pynpact import main


def test_process(patcher, gbkfile, executor):

    mock = patcher.patch('pynpact.main._process')

    assert main.process('extract', gbkfile, executor=executor)
    assert mock.called


def test__process(executor, gbkconfig, patcher):
    "Make sure the _process function handles the generator and calls the executor"
    a = []
    def planner(config):
        a.append((yield (str, 'asdf', [])))
    main._process(planner, gbkconfig, executor)
    assert a[0] == 'asdf'
    assert executor.tasks['asdf'] == ''
