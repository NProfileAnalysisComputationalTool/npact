#!/usr/bin/env python
import glob
import os
import sys
import subprocess
import shutil

pwd = os.path.abspath(os.path.dirname(__file__))
vedir = os.path.abspath(os.path.join(pwd,"ve"))

if os.path.exists(vedir):
    shutil.rmtree(vedir)

virtualenv_support_dir = os.path.abspath(os.path.join(pwd, "requirements", "virtualenv_support"))

ret = subprocess.call(["python", "virtualenv.py", 
                       "--extra-search-dir=%s" % virtualenv_support_dir,
                       "--no-site-packages",
                       "-p", "python2.6",
                       "--never-download",
                       vedir])
if ret: exit(ret)

ret = subprocess.call([os.path.join(vedir, 'bin', 'pip'), "install",
                       "-E", vedir,
                       "--index-url=''",
                       "--requirement",os.path.join(pwd,"requirements/libs.txt")])
if ret: exit(ret)


the_eggs = [os.path.basename(path) for path in
            glob.glob(os.path.join(pwd, "requirements", "eggs", "*.egg"))]
if the_eggs:
    # only try to easy install eggs if there actually are some
    cmd = ([os.path.join(vedir,"bin/easy_install"),
            '-f', os.path.join(pwd,"requirements/eggs/")] +
           the_eggs)
    ret = subprocess.call(cmd)
    exit(ret)
