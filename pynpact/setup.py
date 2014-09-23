#!/usr/bin/env python
import os, os.path, glob

from setuptools import setup, Extension

setup(name='pynpact',
      version='0.4',
      description='Python N-Profile Analysis Computation Tool',
      author='Luciano Brocchieri and Nathan Bird',
      author_email='nathan@acceleration.net',
      url='http://genome.ufl.edu/npact/',
      packages=['pynpact'],
      package_data={
          'pynpact': ['data/*', 'bin/*']
      },
      install_requires=["biopython>=1.57"],
      tests_require=["pytest>=2.4", "mock>=1", "pytest-mock"]
    )
