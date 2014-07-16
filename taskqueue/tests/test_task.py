import pytest

from taskqueue import NoSuchTaskError
from taskqueue.task import Task, async_wrapper
import cPickle

import glob
import sys

from mock import MagicMock


@pytest.fixture()
def anontask():
    return Task(func=str)


@pytest.fixture()
def idtask():
    return Task(func=str, tid="Asdf")


def test_Task(anontask):
    assert repr(anontask)
    anontask.tid = "Asdf"
    assert repr(anontask)
    assert '' is anontask.run()


def test_task_picklability(tmpdir, idtask):
    fooname = tmpdir.join('foo')
    with fooname.open('wb') as tf:
        cPickle.dump(idtask, tf)
    with fooname.open('rb') as tf:
        task2 = cPickle.load(tf)
    assert '' == task2.run()


def test_repeated_pickling(tmpdir):
    "Pickling a task should be a noop when the file exists"
    tid = "Asdf"
    Task(str, tid=tid).pickle(str(tmpdir))
    Task(str, tid=tid).pickle(str(tmpdir))
    assert 1 == len(tmpdir.listdir())


def test_randomid_pickling(tmpdir, anontask):
    anontask.pickle(str(tmpdir))
    assert 1 == len(tmpdir.listdir())
    assert anontask.tid


def test_randomid_pickling_collision(tmpdir, anontask, patcher):
    patcher.patch('taskqueue.task.randomid',
                  MagicMock(side_effect=['asdf', 'fdas']))
    tmpdir.join('asdf').open('w')
    tid = anontask.pickle(str(tmpdir))
    assert 'fdas' == tid
    assert 'fdas' == anontask.tid


def test_task_pickle_unpickle_anon(tmpdir, anontask):
    tid = anontask.pickle(str(tmpdir))
    tsk = Task.unpickle(tid, str(tmpdir))
    assert 1 == len(tmpdir.listdir())
    assert tid == tsk.tid
    assert '' == tsk.run()


def test_task_pickle_unpickle_id(tmpdir, idtask):
    tid = idtask.pickle(str(tmpdir))
    tsk = Task.unpickle(tid, str(tmpdir))
    assert 1 == len(tmpdir.listdir())
    assert tid == tsk.tid
    assert '' == tsk.run()


def test_unpickle_missing(tmpdir):
    with pytest.raises(NoSuchTaskError):
        Task.unpickle('asdf', str(tmpdir))


def test_pickle_with_badid(tmpdir, anontask):
    anontask.tid = "/tmp/foobar"
    anontask.pickle(str(tmpdir))
    assert 1 == len(tmpdir.listdir())
    assert tmpdir.listdir()[0].isfile()
    utsk = Task.unpickle(anontask.tid, str(tmpdir))
    assert utsk.tid == anontask.tid


def test_async_wrapper(tmpdir, idtask):
    stderr = sys.stderr
    assert '' == async_wrapper(idtask, str(tmpdir))
    assert 1 == len(glob.glob(str(tmpdir.join('*.stderr'))))
    assert stderr is sys.stderr
