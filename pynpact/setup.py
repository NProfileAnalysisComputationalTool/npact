#!/usr/bin/env python
from setuptools import setup

setup(name='pynpact',
      version='3.0.0rc1',
      description='Python N-Profile Analysis Computation Tool',
      author='Luciano Brocchieri, Nathan Bird, Russ Tyndall',
      author_email='russ@acceleration.net',
      url='http://genome.ufl.edu/npact/',
      packages=['pynpact'],
      package_data={
          'pynpact': ['data/*', 'bin/*']
      },
      scripts=['bin/pynpact'],
      tests_require=["pytest>=2.4", "mock>=1", "pytest-mock"]
)
