import os.path
import pytest
import py
from pynpact.main import _process

@pytest.fixture()
def plan_processor(executor):
    def func(planner, config):
        return _process(planner, config, executor)
    return func
