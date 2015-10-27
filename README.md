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

# Download

To get the latest stable released version of the code (i.e. what's
running on http://genome.ufl.edu/npact/) go to
"https://github.com/NProfileAnalysisComputationalTool/npact/releases".

to get the latest development version of the code click the "Download
ZIP" button on the right.

# Using

This project is primarily meant to be used as a website (See the
npactflask component). The pynpact component can be used directly on the
command line though this has not been tested recently.

## Requirements

The required libraries are packaged with this
project as much as possible to make it easier to get started and more resistant to
external changes.

You will need [Python](http://python.org/) and [Make](http://www.gnu.org/s/make/)
on your system before getting started.
See [requirements.md](/requirements.md) for more information.


## Building

1. run `bootstrap.py` to get everything going.

        python bootstrap.py

If that completes successfuly then you should see "Successfully
finished bootstrap.py". If it doesn't, the rest won't work.

## Running in Development

### Development mode

There is a development webserver bundled in that will help for local development. To use:

1. `source ve/bin/activate[.csh]`
2. `npactflask/bin/devserver`
3. open http://127.0.0.1:5000/npact/

NB: Apache normally runs as a different user and you may encounter
permissions issues if you run the development server in the same
directory that has been served under Apache.

To log in to the npact management page, the user is `npactmanager`

The password should be located in `.htpasswd` file in the webroot e.g. from the
apache.conf:

    authuserfile /var/www/html/genome.ufl.edu/npact/.htpasswd

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
[Flask](http://flask.pocoo.org). See [npactflask/README.md](/npactflask/README.md).
