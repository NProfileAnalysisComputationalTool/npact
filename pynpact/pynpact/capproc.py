"""This module contains our extensions to subprocess

The name captured proc originally was because a big purpose was to
capture the output and log it. Now it does a whole bunch more than
just log the output.

Be warned, if you see your pipelines hanging read
http://old.nabble.com/subprocess.Popen-pipeline-bug--td16026600.html
and set ``close_fds=True``  (now the default if you haven't set it)

"""
from __future__ import absolute_import
import errno
import logging
import os.path
import subprocess
import threading
import time
from contextlib import contextmanager
from subprocess import PIPE



def check(cmd, **kwargs):
    """Check if a subprocess call exits cleanly.

    Turns off all logging and suppresses all output by default.
    """
    kwargs.setdefault('stdout', False)
    kwargs.setdefault('stderr', False)
    kwargs.setdefault('logger', False)
    return 0 == capturedCall(cmd, **kwargs)


def capturedCall(cmd, check=False, **kwargs):
    """Call capturedPopen and wait for to exit.

    See :py:func:`capturedPopen` for more documentation.

    Additional Args:

    * check (default False): If true then it will check the return
      code of the process and raise a subprocess.CalledProcessError if
      it is non-zero.

    """
    p = capturedPopen(cmd, **kwargs)

    def wait(monitor):
        monitor = getattr(p, 'std%s_monitor_thread' % monitor, None)
        if monitor:
            monitor.join()
    wait('out')
    wait('err')

    rc = p.wait()

    if check and rc != 0:
        # So the CalledProcessError has a bug such that it can't be
        # UNpickled.  http://bugs.python.org/issue1692335#msg133581 by
        # setting the 'args' attribute we workaround that and allow
        # this code to be used inside a multiprocessing.Pool
        args = (rc, p.cmdname)
        cpe = subprocess.CalledProcessError(*args)
        cpe.args = args
        raise cpe
    return rc


def selfCaptured(klass):
    def mungekw(self, kwargs):
        kwargs.setdefault('logger', self.logger)
        return kwargs

    def add(func):
        def newfunc(self, cmd, **kwargs):
            return func(cmd, **mungekw(self, kwargs))
        setattr(klass, func.__name__, newfunc)

    add(capturedCall)
    add(capturedPopen)
    add(guardPopen)
    return klass


class CapProc(object):
    """A superclass to provide some of the captured functions as methods.

    Default some parameters (`logger`) based on the classes property.

    """
    def mungekw(self, kwargs):
        log = (getattr(self, 'logger', None) or
               getattr(self, 'log', None) or
               logging.root)
        kwargs.setdefault('logger', log)
        return kwargs

    def capturedCall(self, cmd, **kwargs):
        return capturedCall(cmd, **self.mungekw(kwargs))

    def capturedPopen(self, cmd, **kwargs):
        return capturedPopen(cmd, **self.mungekw(kwargs))

    def guardPopen(self, cmd, **kwargs):
        return guardPopen(cmd, **self.mungekw(kwargs))

    def guarded_stdout_lines(self, cmd, **kwargs):
        return guarded_stdout_lines(cmd, **self.mungekw(kwargs))


def capturedPopen(cmd, stdin=None, stdout=None, stderr=None,
                  logger=logging.root, cd=None,
                  stdout_level=logging.INFO,
                  stderr_level=logging.WARNING,
                  filter=None,
                  log_command=True,
                  **kwargs):
    """A wrapper around subprocess.Popen that offers the following extensions:

     * stdin, stdout, stderr can be specified as False which will then
       pass an open fd to /dev/null (using os module)

     * if `logger` is provided (default=root logger) then log the
       output of stdout and stderr (only if those streams aren't being
       piped). These will be logged to two new loggers that are
       children of the passed in logger

     * Adds attribute `cmdname` to the returned popen object which is,
       as best as we can divine, the name of the binary being run.

    """
    # We use None as sigil values for stdin,stdout,stderr above so we
    # can distinguish from the caller passing in Pipe.

    if os.name == 'posix' and 'close_fds' not in kwargs:
        # http://old.nabble.com/subprocess.Popen-pipeline-bug--td16026600.html
        kwargs['close_fds'] = True

    if cd:
        # subprocess does this already with the 'cwd' arg, convert cd
        # over so as not to break anyone's code.
        kwargs['cwd'] = cd

    if not isinstance(cmd, basestring):
        cmd = [str(e) for e in cmd]

    if logger and log_command:
        # if we are logging, record the command we're running,
        if 'cwd' in kwargs:
            cwd = " in " + kwargs.get('cwd')
        else:
            cwd = ''
        logger.debug("Running cmd: `%s`%s",
                     (cmd if isinstance(cmd, basestring)
                      else subprocess.list2cmdline(cmd)),
                     cwd)

    # A list of FDs that were opened in the parent, to be passed to
    # child, that need to be closed once that process has been spawned
    close_in_parent = []

    if stdin is False:
        stdin = os.open(os.devnull, os.O_RDONLY)
        close_in_parent.append(stdin)

    def out(arg):
        # figure out what to pass to the stdout stream
        if arg:
            return arg  # specified: use that
        elif arg is False:
            fd = os.open(os.devnull, os.O_WRONLY)
            close_in_parent.append(fd)
            return fd
        elif logger:
            return PIPE
        else:
            return None

    p = subprocess.Popen(cmd, stdin=stdin,
                         stdout=out(stdout),
                         stderr=out(stderr),
                         **kwargs)

    for fd in close_in_parent:
        os.close(fd)

    # try to get a simple name for the command.
    if kwargs.get('shell'):
        p.cmdname = 'sh'
    else:
        p.cmdname = os.path.basename(cmd[0])

    if logger:
        cmdname = p.cmdname[:-4] if p.cmdname.endswith('.exe') else p.cmdname

        def monitor(level, src, name):
            lname = "%s.%s" % (cmdname, name)
            sublog = logger.getChild(lname)
            def tfn():
                l = src.readline()
                while l != "":  # The EOF sigil
                    sublog.log(level, l.rstrip())
                    l = src.readline()

            th = threading.Thread(target=tfn, name=lname)
            p.__setattr__("std%s_monitor_thread" % name, th)
            th.start()

        if stdout is None: monitor(stdout_level, p.stdout, "out")
        if stderr is None: monitor(stderr_level, p.stderr, "err")
    return p


@contextmanager
def guardPopen(cmd, **kwargs):
    """The with-block combination of capturedPopen and guardPopened.

    Accepts in kwargs:

    * keys as defined by :py:func:`capturedPopen`
    * keys as defined by :py:func:`ensure_popen_exits`

    """
    guardargs = {'logger': kwargs.get('logger')}

    def popkey(key):
        if key in kwargs:
            guardargs[key] = kwargs.pop(key)
    popkey('check')
    popkey('timeout')
    popkey('timeout_count')
    with guardPopened(capturedPopen(cmd, **kwargs), **guardargs) as p:
        yield p


@contextmanager
def guardPopened(popen, logger=None, **kwargs):
    """Supervise the given popen process ensuring exist/termination.

    This is a context manager function that will:
     * terminate the process if the exit is due to an exception
       (leaving the exception active).
     * Call :py:func:`ensure_popen_exits` (see f)
    """
    try:
        yield popen
    except Exception:
        terminate_process(popen, logger)
        popen = None  # skip the finally block's cleanup.
        raise
    finally:
        if popen:
            ensure_popen_exits(popen, logger=logger, **kwargs)


def ensure_popen_exits(popen, check=True, timeout=0.4, timeout_count=2,
                       logger=None):
    """Ensure a popen closes one way or another.

    * Wait `timeout` seconds `timeout_count` times for the process to exit
    * terminate the process if it hasn't
    * if check: (default true) raise an exception if the process has
      a non-zero returncode.

    """
    if popen:
        # try waiting a little for the process to finish
        while popen.poll() is None and timeout > 0 and timeout_count > 0:
            if logger:
                cmdname = getattr(popen, 'cmdname', '<unknown>')
                logger.debug(
                    "%s hasn't exited, waiting for %ss (%s tries left)",
                    cmdname, timeout, timeout_count)
            time.sleep(timeout)
            # wait longer each iteration
            timeout = timeout * 2
            timeout_count -= 1

        terminate_process(popen, logger)

        # Did it exit abnormally?
        if check and popen.returncode != 0:
            cmdname = getattr(popen, 'cmdname', '<unknown>')
            raise subprocess.CalledProcessError(popen.returncode, cmdname)


def terminate_process(popen, logger=None, msglevel=logging.FATAL):
    """Do our best to terminate a process checking for windows shenanigans.

    If the process has already exited then do nothing
    """
    # If it hasn't exited yet, assume it's hung.
    if popen and popen.poll() is None:
        if logger and msglevel:
            logger.log(msglevel, "Terminating %s",
                       getattr(popen, 'cmdname', '<unknown>'))
        try:
            popen.terminate()
        except OSError, ose:
            if ose.errno != errno.ESRCH:  # no such process
                raise
        # handle race condition; if the process has exited since the
        # last `poll()`, right now we only see this on
        # windows. Determined error specifics by testing via ipython
        # on windows, trying to `terminate` a completed process.
        except Exception, we:
            if hasattr(we, 'winerror'):
                # 5 is 'Access Denied', could be either our
                # process is dead or the PID is a different
                # process now.
                if we.winerror != 5:
                    raise
                elif logger:
                    logger.fatal(
                        "Windows errno 5 implies proc already terminated: %s",
                        we)
            else:
                # not a windows error, handle normally
                raise


def guarded_stdout_lines(cmd, **kwargs):
    "returns an iterator for stdout lines"
    kwargs['stdout'] = PIPE
    with guardPopen(cmd, **kwargs) as proc:
        try:
            for l in iter(proc.stdout.readline, b''):
                yield l.rstrip()
            proc.wait()
        except GeneratorExit:
            # raised by user of this generator calling .close() (or
            # the generator getting GCed)
            # http://docs.python.org/2/reference/expressions.html#generator-iterator-methods

            # ensure the process has exited.
            terminate_process(
                proc, logger=kwargs.get('logger'), msglevel=False)

            # make proc appear to have exitted normally to bypass any
            # other cleanup.
            proc.returncode = 0
            raise


# TODO: I'm not sure if this is the right abstraction; a class you
# could either call with all the cmds or add to one at a time would be
# nice.
def pipeline(*cmds, **kwargs):
    """Pipe a series of subprocess's output together.

    :param stdin: if given will be the stdin of the first process
    :param stdout: if given will be the stdout of the last process

    Every cmd will be called with capturedPopen, any remaining kwargs
    will be given to every call.

    Return value is the list of popen objects

    """
    if len(cmds) == 0:
        return None
    elif len(cmds) == 1:
        return [capturedPopen(cmds[0], **kwargs)]

    first_stdin = kwargs.pop('stdin', None)
    final_stdout = kwargs.pop('stdout', None)

    popens = []
    for cmd in cmds:
        if cmd is cmds[0]:
            stdin = first_stdin
            stdout = PIPE
        elif cmd is cmds[-1]:
            stdin = popens[-1].stdout
            stdout = final_stdout
        else:
            stdin = popens[-1].stdout
            stdout = PIPE
        popens.append(capturedPopen(cmd, stdin=stdin, stdout=stdout, **kwargs))
    return popens


@contextmanager
def environment(**kwargs):
    "Add extra environment variables as a context manager"
    old_env = {}
    for key, val in kwargs.items():
        old_env[key] = os.environ.pop(key, None)
        if val:
            os.environ[key] = val
    try:
        yield
    finally:
        for key, val in old_env.items():
            if val is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = val

# Copright (c) 2011,2014  Accelerated Data Works
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:

#     Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.

#     Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials provided
#     with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
