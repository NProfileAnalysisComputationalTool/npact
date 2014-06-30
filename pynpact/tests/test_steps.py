import pytest

from pynpact import steps


def test_Task():
    def sum_(a, b):
        return a + b
    assert 3 == steps.Task(sum_, 1, 2)()
    assert 5 == steps.Task(sum_, b=2, a=3)()
    assert type(repr(steps.Task)) == str


def test_delay():
    def sum_(a,b):
        return a+b

    assert hasattr(steps.delay(sum_)(1,2), '__call__')
    assert 3 == steps.delay(sum_)(1,2)()
    assert 3 == steps.delay(sum_)(a=1, b=2)()


@pytest.fixture()
def bs(tmpdir):
    d = {'basename': 'foobar.gbk'}
    e = object()
    return steps.BaseStep(d, str(tmpdir), e)


def test_BaseStep_from(bs):
    bs2 = steps.BaseStep.fromstep(bs)
    assert bs is not bs2
    assert bs.config is bs2.config
    assert bs.outputdir is bs2.outputdir
    assert bs.executor is bs2.executor


def test_BaseStep_derive(bs, tmpdir):
    result = bs.derive_filename('asdf', 'foo')
    assert str(tmpdir.join('foobar-asdf.foo')) == result
