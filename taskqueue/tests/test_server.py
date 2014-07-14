import pytest
from taskqueue import server, NoSuchTaskError
from taskqueue.server import Task, Server, async_wrapper
import cPickle
import tempfile
import glob
import sys

from mock import patch, MagicMock


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


def test_task_picklability(tmpdir, anontask):
    fooname = tmpdir.join('foo')
    with fooname.open('wb') as tf:
        cPickle.dump(anontask, tf)
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
    patcher.patch('taskqueue.server.randomid',
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


def test_server_workdir(patcher, tmpdir):
    patcher.patch('taskqueue.BASE_DIR', str(tmpdir))
    s = Server()
    target = tmpdir / 'queue'
    assert target.exists()
    assert str(target) == str(s.work_dir)


@pytest.fixture()
def aserver(tmpdir):
    return Server(work_dir=str(tmpdir))


def test_instantiation(aserver):
    "Just make sure we're working at all"
    assert aserver


def test_async_wrapper(tmpdir, idtask):
    stderr = sys.stderr
    assert '' == async_wrapper(idtask, str(tmpdir))
    assert 1 == len(glob.glob(str(tmpdir.join('*.stderr'))))
    assert stderr is sys.stderr


def test__enqueue(aserver, idtask):
    async_res = aserver._enqueue(idtask)
    assert aserver.tasks[idtask.tid]
    assert '' == async_res.get()
    assert 1 == len(aserver.tasks)


def test_already__enqueue(aserver, idtask):
    o = object()
    aserver.tasks[idtask.tid] = o
    p = aserver._enqueue(idtask)
    assert o is p
    assert 1 == len(aserver.tasks)


def test_enqueue(aserver, tmpdir):
    tid = aserver.enqueue(str, tid='asdf')
    assert 'asdf' is tid
    assert 1 == len(tmpdir.listdir())
    assert '' is aserver.result(tid)


def test_enqueue_after_missing(aserver):
    with pytest.raises(NoSuchTaskError):
        aserver.enqueue(str, after=['asdf'])
