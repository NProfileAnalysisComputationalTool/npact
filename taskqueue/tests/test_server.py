import pytest
from taskqueue.server import Task, Server, async_wrapper
import cPickle
import tempfile
import glob
import sys


def test_Task():
    def sum_(a, b):
        return a + b
    assert 3 == Task(sum_, [1, 2]).run()
    assert 5 == Task(sum_, [], {'b': 2, 'a': 3}).run()
    assert type(repr(Task)) == str


def test_pickler(tmpdir):
    task = Task(sum, [[1, 2, 3, 4]])
    fooname = tmpdir.join('foo')
    with fooname.open('wb') as tf:
        cPickle.dump(task, tf)
    with fooname.open('rb') as tf:
        task2 = cPickle.load(tf)
    assert 10 == task2.run()


def test_task_properties(tmpdir):
    task_path = tmpdir.join('foo.todo')
    t = Task(sum, [[1, 2]])
    t.path = str(task_path)
    assert str(tmpdir.join('foo')) == t.task_base
    assert False == t.todop
    task_path.write("1")
    assert True == t.todop

    task_path = tmpdir.join('foo')
    t = Task(sum, [[1, 2]])
    t.path = str(task_path)
    assert str(tmpdir.join('foo')) == t.task_base
    assert False == t.todop
    task_path.write("1")
    assert False == t.todop


@pytest.fixture()
def aserver(tmpdir):
    return Server(work_dir=str(tmpdir))


def test_instantiation(aserver):
    "Just make sure we're working at all"
    assert aserver


def test_pickle_task(aserver):
    tid, path = aserver.pickle_task(Task(sum))
    assert tid
    assert path

    t2 = aserver.unpickle_task(tid)
    assert sum == t2.func
    assert tid == t2.tid
    assert path == t2.path


@pytest.fixture()
def atask(tmpdir):
    t = Task(sum, [[1, 2]])
    t.tid = "35415"
    t.path = str(tmpdir.join(t.tid + ".todo"))
    return t


def test_mark_finished(tmpdir, atask):
    atask.mark_finished()
    assert False == atask.todop
    assert 0 == len(glob.glob(str(tmpdir.join('*.todo'))))


def test_async_wrapper(tmpdir, atask):
    assert 3 == async_wrapper(atask)
    assert 1 == len(glob.glob(str(tmpdir.join('*.stderr'))))
    assert False == atask.todop


def test__enqueue(aserver, atask):
    async_res = aserver._enqueue(atask)
    assert aserver.tasks[atask.tid]
    assert 3 == async_res.get()
