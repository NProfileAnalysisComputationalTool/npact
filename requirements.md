### External requirements

These will need to be setup on your system before getting started.

* [Python](http://python.org/): The bulk of the glue code is written
  in Python. It is targeting python 2.6 though 2.7 should work just as
  well. Python 3 compatibility has not been tested - and probably won't
  work due to biopython.
* A [posix](http://en.wikipedia.org/wiki/POSIX) environment for
  python. It works on CentOS and Ubuntu; it should work on any posix
  environment build of python. Mac OS X is expected to work.
* A C compiler: the actual analysis code is written in C. Tested with
  gcc 4.8.4. Others should work, - I don't think there is anything too
  crazy being used.
* [Make](http://www.gnu.org/s/make/): A makefile is used to build all
  the C.
* A PostScript viewer to view the output files.
* [Git](http://git-scm.com/) (OPTIONAL): The version control system
  this project is maintained in. Will be necessary to record changes
  but not for running the project.
* [ps2pdf](http://ghostscript.com/doc/current/Ps2pdf.htm) (OPTIONAL):
  If present NPACT will convert the rendered Postscript to a PDF file
  using this program.

### Packaged requirements

The system uses several python packages for deployment:

* [virtualenv](http://pypi.python.org/pypi/virtualenv)
* [distribute](http://pypi.python.org/pypi/distribute)
* [pip](http://www.pip-installer.org/en/latest/index.html)

These are already included.

The `requirements.txt` contains the exact libraries beyond that. It is
in a [format][req-file-format] that pip understands.

Notes:

* Biopython (1.60): Used to query entrez and read information out of
  GenBank files. This can probably be upgraded without any hassle.
* Flask (0.10): Used to build the website interface.

[req-file-format]: https://pip.pypa.io/en/latest/reference/pip_install.html#requirements-file-format
