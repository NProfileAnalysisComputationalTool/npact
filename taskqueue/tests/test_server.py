import pytest
import taskqueue
from taskqueue import server, NoSuchTaskError
from taskqueue.server import Task, Server, async_wrapper
import cPickle
import tempfile
import glob
import sys
import logging


def test_server_workdir(patcher, tmpdir):
    patcher.patch('taskqueue.BASE_DIR', str(tmpdir))
    s = Server()
    target = tmpdir / 'queue'
    assert target.exists()
    assert str(target) == str(s.work_dir)


@pytest.fixture()
def idtask(magmock):
    magmock.tid = 'Asdf'
    return magmock


@pytest.fixture()
def aserver(tmpdir, magmock):
    return Server(work_dir=str(tmpdir), pool=magmock)


def test_instantiation(aserver):
    "Just make sure we're working at all"
    assert aserver


def test__enqueue(aserver, idtask):
    aserver._enqueue(idtask)
    assert aserver.tasks[idtask.tid]
    assert aserver.pool.apply_async.called
    assert aserver.get_task(idtask.tid)


def test_enqueue(aserver, tmpdir, patcher):
    patcher.patch_object(aserver, '_enqueue')
    tid = aserver.enqueue(str, tid='asdf')
    assert 'asdf' is tid
    assert 1 == len(tmpdir.listdir())
    assert aserver._enqueue.called


def test_enqueued_already(aserver, patcher):
    patcher.patch_object(aserver, '_enqueue')
    aserver.tasks['asdf'] = object()
    assert 'asdf' == aserver.enqueue(str, tid='asdf')
    assert not aserver._enqueue.called


def test_enqueue_after_missing(aserver):
    with pytest.raises(NoSuchTaskError):
        aserver.enqueue(str, after=['asdf'])


@pytest.fixture()
def aRunningServer(request, patcher):
    patcher.patch('taskqueue.client.ENSURE_DAEMON', False)

    sm = taskqueue.get_ServerManager(make_server=True)
    logging.info("Opening a socket at %s", sm.address)
    sm.start()
    request.addfinalizer(sm.shutdown)
    return sm.Server()


def test_enqueue_after(aRunningServer, tmpdir):
    tid1 = aRunningServer.enqueue(str)
    tid2 = aRunningServer.enqueue(str, after=[tid1])
    assert '' == aRunningServer.get_task(tid2).get(1)
