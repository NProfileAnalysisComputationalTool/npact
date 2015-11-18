#!/usr/bin/env python
from setuptools import setup

setup(name='pynpact',
      version='0.6',
      description='Python N-Profile Analysis Computation Tool',
      author='Luciano Brocchieri and Nathan Bird',
      author_email='nathan@acceleration.net',
      url='http://genome.ufl.edu/npact/',
      packages=['pynpact'],
      package_data={
          'pynpact': ['data/*', 'bin/*']
      },
      scripts=['bin/pynpact'],
      install_requires=["biopython>=1.57"],
      tests_require=["pytest>=2.4", "mock>=1", "pytest-mock"]
)
