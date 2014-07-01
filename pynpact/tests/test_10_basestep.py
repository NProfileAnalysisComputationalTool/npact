import pytest

from pynpact import steps


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
