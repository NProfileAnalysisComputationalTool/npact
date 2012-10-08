#!/usr/bin/env python
from distutils.core import setup
import os, os.path


def recursive_files(*roots):
    for r in roots:
        for root,directories,files in os.walk(r):
            for i in files:
                yield os.path.join(root,i)

setup(name='npactweb',
      version='0.4',
      description='Website for PYNPACT, the Python N-Profile Analysis Computation Tool',
      author='Nathan Bird',
      author_email='nathan@acceleration.net',
      url='http://genome.ufl.edu/npact/',
      py_modules=['settings'],
      packages=['npactweb'],
      package_data={'npactweb': list(recursive_files('static','templates'))},
      requires=["biopython(>=1.57)",
                "pynpact",
                "django(==1.3)"],
      scripts=['cleanup.py', 'django-main.fcgi', 'django-process.fcgi']
     )
