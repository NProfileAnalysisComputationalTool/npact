import os, os.path, logging, subprocess, threading, time, errno, hashlib
import tempfile
from subprocess import PIPE
from contextlib import contextmanager



def reducehashdict(dict,keys) :
    """pull the given keys out of the dictionary, return the reduced
    dictionary and the sha1 hash of that set of key values.
"""
    outdict= {}
    h = hashlib.sha1()
    for k in sorted(keys) :
        val = dict.get(k)
        if val is not None:
            h.update(k)
            h.update(str(val))
            outdict[k]=val

    if len(outdict) :
        return outdict,h.hexdigest()
    else :
        return {},None

def ensure_dir(dir, logger=None) :
    if not os.path.exists(dir) :
        try:
            if logger: logger.debug("Making dir: %s", dir)
            os.makedirs(dir)
            if logger: logger.info("Created dir: %s", dir)
        except OSError, e :
            #not entirely sure why we are getting errors,
            #http://docs.python.org/library/os.path.html#os.path.exists
            #says this could be related to not being able to call
            #os.stat, but I can do so from the command line python
            #just fine.
            if os.path.exists(dir) :
                if logger: logger.debug("Erred, already exists: e.errno: %s",e.errno)
                return
            else:
                raise


def withDir(dir, fn, *args,**kwargs) :
    olddir = os.getcwd()
    try :
        os.chdir(dir)
        return fn(*args,**kwargs)
    finally :
        os.chdir(olddir)


def pprint_bytes(bytes) :
    suffix = 'B'
    bytes= float(bytes)
    if bytes >= 1024 :
        bytes = bytes / 1024
        suffix = 'KB'
    if bytes >= 1024 :
        bytes = bytes / 1024
        suffix = 'MB'
    if bytes >= 1024 :
        bytes = bytes / 1024
        suffix = 'GB'
    if bytes >= 1024 :
        bytes = bytes / 1024
        suffix = 'TB'
    return '%.2f%s' % (bytes,suffix)




def exec_external( proc ):
    return subprocess.call(['python', proc])


def capturedCall(cmd, **kwargs) :
    """Do the equivelent of the subprocess.call except
    log the stderr and stdout where appropriate."""
    p= capturedPopen(cmd,**kwargs)
    rc = p.wait()
    #this is a cheap attempt to make sure the monitors
    #are scheduled and hopefully finished.
    time.sleep(0.01)
    time.sleep(0.01)
    return rc


def capturedPopen(cmd, stdin=None, stdout=None, stderr=None,
                  logger=logging,
                  stdout_level=logging.INFO,
                  stderr_level=logging.WARNING, **kwargs) :
    """Equivalent to subprocess.Popen except log stdout and stderr
    where appropriate. Also log the command being called."""
    #we use None as sigil values for stdin,stdout,stderr above so we
    # can distinguish from the caller passing in Pipe.

    if not kwargs.has_key("close_fds") and os.name == 'posix' :
        #http://old.nabble.com/subprocess.Popen-pipeline-bug--td16026600.html
        kwargs['close_fds'] = True

    if not isinstance(cmd,str) :
        cmd = [str(e) for e in cmd]

    if(logger):
        #if we are logging, record the command we're running,
        #trying to strip out passwords.
        logger.debug("Running cmd: %s",
                     isinstance(cmd,str) and cmd or subprocess.list2cmdline(cmd))

    p = subprocess.Popen(cmd, stdin=stdin,
                         stdout=(stdout or (logger and PIPE)),
                         stderr=(stderr or (logger and PIPE)),
                         **kwargs)
    if logger :
        def monitor(level, src, name) :
            #if the cmd[0] (the binary) contains a full path, just get the name
            lname = "%s.%s" % (os.path.basename(cmd[0]), name)
            if(hasattr(logger, 'name')) :
                lname = "%s.%s" % (logger.name, lname)
            sublog = logging.getLogger(lname)

            def tfn() :
                l = src.readline()
                while l != "":
                    sublog.log(level,l.strip())
                    l = src.readline()

            th = threading.Thread(target=tfn,name=lname)
            p.__setattr__("std%s_thread" % name, th)
            th.start()

        if stdout == None : monitor(stdout_level, p.stdout,"out")
        if stderr == None : monitor(stderr_level, p.stderr,"err")
    return p


@contextmanager
def guardPopen(cmd, timeout=0.1, timeout_count=2, **kwargs) :
    popen = None
    logger = kwargs.get('logger')
    try :
        popen = capturedPopen(cmd,**kwargs)
        yield popen
    finally :
        if popen:
            while popen.poll() == None and timeout and timeout_count > 0 :
                if logger:
                    logger.debug("Things don't look right yet waiting for %ss (%s tries left)", 
                                 timeout, timeout_count)
                time.sleep(timeout)
                timeout_count -=1
            if popen.poll() == None:
                if logger:
                    logger.exception("Terminating %s", cmd[0])
                popen.terminate()

def selfCaptured(klass) :
    def mungekw(self,kwargs) :
        if(not kwargs.has_key("logger")): kwargs["logger"] = self.logger
        return kwargs
    def add(func) :
        def newfunc(self,cmd,**kwargs) :
            return func(cmd, **mungekw(self,kwargs))
        setattr(klass,func.__name__,newfunc)

    add(capturedCall)
    add(capturedPopen)
    add(guardPopen)
    return klass

def which(program):
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def stream_to_file(stream,path,bufsize=8192) :
    with open(path, "wb") as f:
        bytes=0
        while True :
            buf = stream.read(bufsize)
            if buf == "" :  break #EOF
            bytes += len(buf)
            f.write(buf)
        return bytes



@contextmanager
def mkstemp_overwrite(destination, conflict_overwrite=True, logger=None,
                      cleanup=True, **kwargs) :
    """For writing to a temporary file and then move it ontop of a
    (possibly) existing file only when finished.  This enables us to
    perform long running operations on a file that other people might
    be using and let everyone else see a consistent version

 * conflict_overwrite=True: whether or not to overwrite the file if
   it's been modified in the intermediate time.
 * logger=None: if provided log information to it.
 * cleanup=True: ensure the tempfile is delete when we exit the
   function (e.g. an error or conflict)
 * other args are passed to tempfile.mkstemp

Example:
with mkstemp_overwrite('foobar.txt') as f:
   #long time intensive processing
   f.write('stuff\n')


    """

    mtime1 = mtime2 = None
    if os.path.exists(destination) :
        mtime1 = os.path.getmtime(destination)

    (fd,path) = tempfile.mkstemp(**kwargs)
    try :
        filelike = os.fdopen(fd,'wb')
        yield filelike
        filelike.close()

        if os.path.exists(destination) :
            mtime2 = os.path.getmtime(destination)

        if mtime1 != mtime2 and logger:
            logger.warning("Potential conflict on %r, overwrite: %s; ts1:%s, ts2:%s",
                               destination, conflict_overwrite, mtime1, mtime2)

        if conflict_overwrite or mtime1 == mtime2 :
            #TODO: permissions?
            os.rename(path,destination)
    finally :
        if cleanup and os.path.exists(path) :
            if logger :
                logger.info("Cleaning up leftover tempfile %r", path)
            os.remove(path)


def is_outofdate(filename, *dependencies) :
    if not os.path.exists(filename) : return True

    mtime = os.path.getmtime(filename)
    return any(os.path.getmtime(d) > mtime for d in dependencies if d)

def derivative_filename(base, part, replace_ext=True, outputdir=None, dependencies=[]) :
    """Build the filename of a derivative product of the original
    file. If the derivative file already exists return whether it
    is out of date"""

    if not part[0] == "." :
        part = "." + part

    if outputdir is None :
        outputdir = os.path.dirname(base)
    filename = os.path.basename(base)

    if replace_ext :
        filename = os.path.splitext(filename)[0]

    outfilename = os.path.join(outputdir, filename + part)

    return outfilename

def safe_produce_new(outfilename, func,
                     replace_ext=True, force=False,dependencies=[], **kwargs) :
    outofdate = is_outofdate(outfilename, *dependencies)
    logger=kwargs.get('logger')
    if outofdate or force:
        if logger:
            logger.debug("Regenerating, checked:%d force:%r", len(dependencies),force)
            
        with mkstemp_overwrite(outfilename,**kwargs) as f :
            func(f)
    elif kwargs.get('logger',False) :
        kwargs.get('logger').debug("Skipped producing %r", outfilename)
    return outfilename

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
