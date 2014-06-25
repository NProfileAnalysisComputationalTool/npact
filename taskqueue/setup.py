#!/usr/bin/env python
from setuptools import setup
from setuptools.command.test import test as TestCommand
import sys


class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = None

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

description = """
Basic process pool daemon that can run python tasks asynchronously
""".strip()

setup(name='taskqueue',
      version='0.2.1',
      description=description,
      author='Nathan Bird',
      author_email='nathan@acceleration.net',
      url='http://genome.ufl.edu/npact/',
      packages=['taskqueue'],
      requires=["lockfile"],
      tests_require=["pytest>=2.4", "mock>=1"],
      scripts=['tqdaemon.py'],
      cmdclass={'test': PyTest})
