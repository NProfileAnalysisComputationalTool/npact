#!/usr/bin/env python

from optparse import OptionParser
import glob
import os
import sys
import subprocess
import shutil
import logging

pwd = os.path.abspath(os.path.dirname(__file__))
vedir = os.path.abspath(os.path.join(pwd,"ve"))

def cleanup_existing():
    if os.path.exists(vedir):
        logging.debug("Attempting to shutdown tqdaemon")
        rc = subprocess.call(['ve/bin/python', 'manage.py', 'tqdaemon', 'stop'])
        if rc != 0:
            logging.info("Unable to shutdown tqdaemon")

        logging.info("Removing existing virtual environment")
        shutil.rmtree(vedir)

def perform_checkouts():
    logging.debug("Checking out submodules")
    ret = subprocess.call(['git', 'submodule', 'update', '--init'])
    if ret:
        print "Error updating git libraries; may not build properly."

def init_virtualenv():
    logging.info("Creating virtualenvironment")
    virtualenv_support_dir = os.path.abspath(os.path.join(pwd, "lib", "virtualenv_support"))
    ret = subprocess.call(["python", "virtualenv.py",
                           "--extra-search-dir=%s" % virtualenv_support_dir,
                           "--distribute",
                           "--no-site-packages",
                           "--prompt=(npact)",
                          "-p", "python2.6",
                          "--never-download",
                           vedir],
                          cwd=pwd)
    if ret:
        logging.critical("Failed creating virtualenv: %s", ret)
        exit(ret)

    logging.debug("Installing libraries")
    ret = subprocess.call([os.path.join(vedir, 'bin', 'pip'), "install",
                           "-E", vedir,
                           "--index-url=''",
                           "--requirement",os.path.join(pwd,"requirements.txt")],
                          cwd=pwd)
    if ret:
        logging.critical("Failed installing libraries from requirements.txt")
        exit(ret)

def create_aux_directories():
    if not os.path.exists('webroot'):
        os.makedirs('webroot')

def build_pynpact():
    logging.info("Building pynpact C code.")
    pynpact_dir = os.path.join(os.path.dirname(__file__), "pynpact/")
    ret = subprocess.call(["make"], cwd=pynpact_dir)
    if ret:
        print "Make failed C programs may not be available."
    else:
        print "Successfully finished bootstrap.py"

if __name__ == '__main__':
    parser = OptionParser("""%prog [options]

    Runs the initialization and building code to get the system into a
    working state after a fresh checkout.
    
    """)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")
    (options,args) = parser.parse_args()
    logging.basicConfig(level=(options.verbose and logging.DEBUG or logging.INFO),
                        format="%(asctime)s %(levelname)-8s %(message)s",
                        datefmt='%H:%M:%S')
    cleanup_existing()
    perform_checkouts()
    init_virtualenv()
    create_aux_directories()
    build_pynpact()
    
