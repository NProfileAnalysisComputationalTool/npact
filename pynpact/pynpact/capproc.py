"""
This module contains our extensions to subprocess
"""
import subprocess
import time
import logging
import os, os.path
from contextlib import contextmanager
from subprocess import PIPE
import threading



def capturedCall(cmd, check=False, **kwargs):
    """
    Call capturedPopen and wait for to exit. 

    See capturedPopen for more documentation.

    Additional Args:
     check (default False):
       If true then it will check the return code of the process and
       raise a subprocess.CalledProcessError if it is non-zero.
    """
    p= capturedPopen(cmd,**kwargs)

    def wait(monitor):
        monitor = getattr(p, 'std%s_monitor_thread' % monitor, None)
        if monitor:
            monitor.join()
    wait('out')
    wait('err')

    rc = p.wait()

    if check and rc != 0:
        raise subprocess.CalledProcessError(rc, p.cmdname)
    return rc


def selfCaptured(klass):
    def mungekw(self,kwargs):
        kwargs.setdefault('logger', self.logger)
        return kwargs
    def add(func):
        def newfunc(self,cmd, **kwargs):
            return func(cmd, **mungekw(self, kwargs))
        setattr(klass, func.__name__, newfunc)

    add(capturedCall)
    add(capturedPopen)
    add(guardPopen)
    return klass



def capturedPopen(cmd, stdin=None, stdout=None, stderr=None,
                  logger=logging,cd=None,
                  stdout_level=logging.INFO,
                  stderr_level=logging.WARNING,
                  filter=None,
                  log_command=True,
                  **kwargs):
    """A wrapper around subprocess.Popen that offers the following extensions:

     * stdin, stdout, stderr can be specified as False which will then
       pass an open fd to /dev/null (using os module)

     * if logger is provided (default=logging) then log the output of
       stdout and stderr (only if those streams aren't being piped). These will be logged to two new loggers
       that are children of the passed in loggerIf a
       'filter' is provided that will be added to the created loggers.

     * Adds attribute cmdname which is, as best as we can divine, the name
       of the binary being run.
   """
    #we use None as sigil values for stdin,stdout,stderr above so we
    # can distinguish from the caller passing in Pipe.

    if os.name == 'posix' and not kwargs.has_key("close_fds"):
        #http://old.nabble.com/subprocess.Popen-pipeline-bug--td16026600.html
        kwargs['close_fds'] = True

    if not isinstance(cmd,str):
        cmd = [str(e) for e in cmd]
        
    if logger and log_command:
        #if we are logging, record the command we're running,
        #trying to strip out passwords.
        logger.debug("Running cmd: %s", 
                     cmd if isinstance(cmd, str) else subprocess.list2cmdline(cmd))
    if cd:
        #subprocess does this already with the 'cwd' arg,
        #convert cd over so as not to break anyone's code.
        kwargs['cwd']=cd

    # A list of FDs that were opened in the parent, to be passed to
    # child, that need to be closed once that process has been spawned
    close_in_parent = []

    if stdin is False:
        stdin=os.open(os.devnull, os.O_RDONLY)
        close_in_parent.append(stdin)

    def out(arg):
        #figure out what to pass to the stdout stream
        if arg            : return arg #specified: use that
        elif arg is False : 
            fd = os.open(os.devnull, os.O_WRONLY)
            close_in_parent.append(fd)
            return fd
        elif logger       : return PIPE
        else              : return None

    p = subprocess.Popen(cmd, stdin=stdin,
                         stdout=out(stdout),
                         stderr=out(stderr),
                         **kwargs)

    for fd in close_in_parent:
        os.close(fd)

    #try to get a simple name for the command.
    if isinstance(cmd, str):
        p.cmdname = os.path.basename(cmd.split(' ')[0])
    else:
        p.cmdname = os.path.basename(cmd[0])

    if logger:
        def monitor(level, src, name):
            lname = "%s.%s" % (p.cmdname, name)
            if logger is logging:
                sublog = logging.getLogger(lname)
            else: 
                sublog = logging.getLogger('%s.%s' % (logger.name,lname))
            if filter:
                try:
                    map(sublog.addFilter, filter)
                except TypeError:
                    if hasattr(filter, 'filter'):
                        sublog.addFilter(filter)

            def tfn():
                l = src.readline()
                while l != "":
                    sublog.log(level, l.strip())
                    l = src.readline()

            th = threading.Thread(target=tfn, name=lname)
            p.__setattr__("std%s_monitor_thread" % name, th)
            th.start()

        if stdout is None: monitor(stdout_level, p.stdout,"out")
        if stderr is None: monitor(stderr_level, p.stderr,"err")
    return p

@contextmanager
def guardPopen(cmd, timeout=0.1, timeout_count=2, **kwargs):
    """Calls capturedPopen and gives the result to guardPopened which
makes sure the process exits as expected.
"""
    with guardPopened(capturedPopen(cmd, **kwargs),
                     timeout=timeout, timeout_count=timeout_count,
                     logger=kwargs.get('logger')) as p:
        yield p


#be warned, if you see your pipelines hanging:
#http://old.nabble.com/subprocess.Popen-pipeline-bug--td16026600.html
#close_fds=True

@contextmanager
def guardPopened(popen, timeout=0.1, timeout_count=2, logger=None):
    """Supervise the given popen process ensuring exist/termination.

    This is a context manager function that will:
     * terminate the process if the exit is due to an exception
       (leaving the exception active).
     * terminate the process and raise an exception if the
       process doesn't exit.
     * raise an exception if the process stops, but has a
       non-zero returncode.
    """
    cmdname = getattr(popen, 'cmdname', '<unknown>')
    try:
        yield popen
    except Exception, e:
        if popen.poll() is None:
            if logger:
                logger.fatal("Terminating %r while handling exception %r)", cmdname, e)
            popen.terminate()
            popen = None  #skip the finally block's cleanup.
        raise
    finally:
        if popen:
            # try waiting a little for the process to finish
            while popen.poll() is None and timeout > 0 and timeout_count > 0:
                if logger:
                    logger.debug("%s didn't exit as expected, waiting for %ss (%s tries left)",
                                 cmdname, timeout, timeout_count)
                time.sleep(timeout)
                timeout_count -= 1

            # If it hasn't exited yet, assume it's hung.
            if popen.poll() is None:
                if logger: logger.fatal("Terminating %s", cmdname)
                popen.terminate()

            # Did it exit abnormally?
            if popen.poll() != 0:
                raise subprocess.CalledProcessError(popen.returncode, cmdname)

# Copright (c) 2011  Accelerated Data Works
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
