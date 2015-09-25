# About

NPACT (N-Profile Analysis and Computation Tool) implements methods to
identify in nucleotide sequences regions of statistically-significant
compositional three-base periodicity associated to ORF structures,
to compare the newly-identified ORFs with pre-annotated collections
of genes and to highlight potential new genes not included in the annotation.
The primary output is a graphical representation of sequence features
represented by frame-specific compositional profiles (by default, GC content)
and sequence segments with significant periodicity (Hits) together with
representation of the pre-annotated genes and, in a separate track, newly
identified ORFs not included in the annotation.

NPACT accepts input sequences of any length and its dynamical graphical
representations allow easy browsing and analysis of an entire prokaryotic genome.

This NPACT project was developed in 2011-2015 by Acceleration.net under
contract from Luciano Brocchieri at the [University of Florida
Genetics Institute](http://www.ufgi.ufl.edu/) based on sequence-analysis
methods and computational tools developed by Luciano Brocchieri at the
[University of Florida Genetics Institute](http://www.ufgi.ufl.edu/).

# Using

This project is primarily meant to be used as a website (See the
npactflask component). The pynpact component can be used directly on the
command line though this has not been tested recently.

## Requirements

As much as possible the required libraries are packaged with this
project to make it easier to get started and more resistant to
external changes.

### External requirements

These will need to be setup on your system before getting started.

* [Python](http://python.org/): The bulk of the glue code is written
  in Python. It is targeting python 2.6 though 2.7 should work just as
  well. Python 3 compatibility has not been tested--and probably won't
  work due to biopython.
* A [posix](http://en.wikipedia.org/wiki/POSIX) environment for
  python. It works on CentOS and Ubuntu; it should work on any posix
  environment build of python. Mac OS X is expected to work.
* A C compiler: the actual analysis code is written in C. Tested with
  gcc 4.8.4. Others should work, I don't think there is anything too
  crazy being used.
* [Make](http://www.gnu.org/s/make/): A makefile is used to build all
  the C.
* Some sort of PostScript viewer to view the output files.
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


## Building

1. run bootstrap.py to get everything going.

        python bootstrap.py

If that completes successfuly then you should see "Successfully
finished bootstrap.py". If it doesn't the rest won't work.

## Running in Development

### Development mode

There is a development webserver bundled in that will help for local development. To use:

1. `source ve/bin/activate[.csh]`
2. `npactflask/bin/devserver`
3. open [http://127.0.0.1:5000/npact/]()

NB: Apache normally runs as a different user and you may encounter
permissions issues if you run the development server in the same
directory that has been served under Apache.

To log in to the npact management page use:
user: npactmanager

The password should be located in .htpasswd file in the webroot
apache.conf:    authuserfile /var/www/html/genome.ufl.edu/npact/.htpasswd

### command line interface

The project was developed with this in mind but has been extended
beyond it so there is no good way to provide detailed configuration
via the command line--it's only available via the website. This is
still useful to run the analysis with all the default options.

1. In a shell from the project root activate the virtualenv:
   `. ve/bin/activate`
2. `pynpact <gbkfile>`
3. It will print a message reporting the output.


# Project Components

In the root of the project are several configuration and build files.

* `requirements.txt` - the list of python packages that are used. Most
  of the referenced packages exist in the `lib` folder; the rest are
  described below.
* `apache.conf` - sample apache configuration (and what's actually
  used on genome.ufl.edu)
* `bootstrap.py` - the script that helps get a new checkout of the
  project running. Run this after a fresh checkout to get the virtual
  environment and dependencies setup.
* `publish.sh` - used for publishing to genome.ufl.edu (or potentially
  another host)
* `lib/` - contains the libraries this project depends on.
* `ve/` - holds the python virtual environment-- won't exist until
  after `bootstrap.py` has been called. Can be deleted and recreated
  with `boostrap.py` at any point.
* `webroot/` - a container for the dynamic portion that will change
  while running as a website. Uploaded files, logs, daemon pidfile,
  and deferred tasks are all written underneath here.

## pynpact/

For Python N-Profile Analysis and Computation Tool.  This is the main
module that interacts with
[Genbank](http://www.ncbi.nlm.nih.gov/genbank/) files to compute the
graphs of N-Profile ratios.

The C code at `pynpact/src` does most of the analysis and the python
code in `pynpact/pynpact` helps glue it together into an easy-to-run
process.

* `pynpact/pynpact/parsing.py` this is used to investigate a GBK file
  and come up with the default options that will be used on a run. To
  change the default analysis options.on the websites run page, go
  here. To change the help text about one of the options, go here.
* `pynpact/pynpact/entrez.py` is used for querying entrez to find gbk
  files.
* `pynpact/pynpact/main.py` actually coordinates the running of the C
  programs and produces the PS output.


## npactflask/

This is the code for the website. It is built on top of
[Flask](http://flask.pocoo.org). See [npactflask/README.md]() in that
folder.


# Contributing

This project uses Git for source control. Git is an excellent system
with lots of good documentation. One of its advantages is that it can
be run disconnected from others and then used to bring disparate
parties back into sync at some later point. See the section "Recording
changes" below for a super quick introduction and links to tutorials.

## Virtualenv

In order for any editing to work well you should activate the python
virtual environment. This will alter the paths that the system looks
for the python and executable code in. From a command line shell in
the root of the project simply type.

    source ve/bin/activate

## Editing the website

Website content is mostly in the templates at
`npactflask/npactflask/templates`. This uses the
[Flask Jinja2 templating ](http://flask.pocoo.org/docs/0.10/templating/). To
alter the default values or help text given on the run page go look in
`pynpact/pynpact/parsing.py`.

This can be tested out by running the development server. See "Using >
Running > Development mode" section above.

## Editing the C code

All of the C lives in `pynpact/src`.

### Building with make

Make is required to build the C code associated with this
project. Make will automatically detect which files have changed since
it was last invoked and only build those pieces.

    cd pynpact
    make


## Recording changes

Frequently recording your changes makes it easier to merge with other
developers later or to track down when a bug or confusing line of code
was introduced.

The simplest form is:

    git commit -m "<A brief descriptive message>"  <the file or files to commit>

How frequently should you do this? Any time you are at a good little
stopping point: e.g.

    git commit -m "Got most of the way through first draft of the README file." README

Later on you can see the entire history of a file:

    git log README

If you are in the middle of editing a file and want to see what all
you've changed?

    git diff <file>

Want to get rid of those changes and go back to the recorded version
that you know worked?

    git checkout <file>

There are many very good tutorials out there;
[Pro Git](http://progit.org/book/) is fantastic and all on
[this list](http://sixrevisions.com/resources/git-tutorials-beginners/)
look good.
