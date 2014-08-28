import logging
import pytest
import mock as mocklib

@pytest.fixture(scope="module", autouse=True)
def setup_logging():
    logging.root.setLevel(logging.DEBUG)
    logging.basicConfig()


class Patcher(object):
    def __init__(self):
        self._patchers = []

    def _stop(self):
        for patcher in self._patchers:
            patcher.stop()

    def _patch(self, method, *args, **kwds):
        if method:
            patch_fn = getattr(mocklib.patch, method)
        else:
            patch_fn = mocklib.patch

        patcher = patch_fn(*args, **kwds)
        self._patchers.append(patcher)

        return patcher.start()

    def patch(self, *args, **kwds):
        return self._patch(None, *args, **kwds)

    def patch_object(self, *args, **kwds):
        return self._patch('object', *args, **kwds)

    def patch_dict(self, *args, **kwds):
        return self._patch('dict', *args, **kwds)

    def patch_multiple(self, *args, **kwds):
        return self._patch('multiple', *args, **kwds)


@pytest.fixture
def patcher(request):
    p = Patcher()
    request.addfinalizer(p._stop)
    return p


@pytest.fixture
def magmock(request):
    return mocklib.MagicMock()
