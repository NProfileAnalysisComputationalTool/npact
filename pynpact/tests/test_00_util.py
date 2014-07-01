from pynpact import util


def test_Task():
    def sum_(a, b):
        return a + b
    assert 3 == util.Task(sum_, 1, 2)()
    assert 5 == util.Task(sum_, b=2, a=3)()
    assert type(repr(util.Task)) == str


def test_delay():
    def sum_(a, b):
        return a+b

    assert hasattr(util.delay(sum_)(1, 2), '__call__')
    assert 3 == util.delay(sum_)(1, 2)()
    assert 3 == util.delay(sum_)(a=1, b=2)()


def test_Hasher():
    h = util.Hasher()
    d = h.hexdigest()
    h2 = h.hash('asdf')
    assert h is h2
    assert d != h.hexdigest()
    d = h.hexdigest()

    h2 = h.hashfiletime(__file__)
    assert h is h2
    assert d != h.hexdigest()
    d = h.hexdigest()

    h2 = h.hashdict({'a': 1, 'b': 2, 'c': None})
    assert h is h2
    assert d != h.hexdigest()
    d = h.hexdigest()

    h2 = h.hashlist([1, 2, 3])
    assert h is h2
    assert d != h.hexdigest()
    d = h.hexdigest()


def test_mkstemp_overwrite(tmpdir):
    filename = tmpdir.join('file')
    with util.mkstemp_overwrite(str(filename)) as out:
        out.write('foo')
        assert len(tmpdir.listdir()) == 1
        assert not filename.exists()
    assert filename.exists()


def test_mkstemp_overwrite2(tmpdir):
    filename = tmpdir.join('file')
    with util.mkstemp_overwrite(str(filename), dir=str(tmpdir)) as out:
        out.write('foo')
        out.flush()
        assert len(tmpdir.listdir()) == 1
        assert not filename.exists()
    assert filename.exists()
