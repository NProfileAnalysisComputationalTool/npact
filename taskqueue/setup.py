#!/usr/bin/env python
from distutils.core import setup
import os, os.path


def recursive_files(*roots):
    for r in roots:
        for root,directories,files in os.walk(r):
            for i in files:
                yield os.path.join(root,i)

setup(name='taskqueue',
      version='0.1',
      description='Basic process pool daemon that can run python tasks asynchronously',
      author='Nathan Bird',
      author_email='nathan@acceleration.net',
      url='http://genome.ufl.edu/npact/',
      packages=['taskqueue'],
      requires=["pynpact",],
      scripts=['tqdaemon.py']
     )
