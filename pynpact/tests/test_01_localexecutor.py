"""A simple executor that will be used for testing.

We need to test this first to make sure if there is an error with the
executor we catch it before it shows up as an error in another test.

"""
import pytest

from pynpact.util import delay


def test_InlineExecutor(executor):
    tid = executor.enqueue(delay(sum)([1, 2, 3]))
    assert 6 == executor.result(tid)

    tid = executor.enqueue(delay(sum)([1, 2, 3]), tid="12414214", after=[tid])
    assert "12414214" == tid
    assert 6 == executor.result(tid)

    with pytest.raises(Exception):
        executor.enqueue(delay(sum)([1]), after="nonexistant")
