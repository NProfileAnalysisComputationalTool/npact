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
