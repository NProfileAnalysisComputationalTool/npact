import time

class Timeout(Exception):
    def __init__(self, **args):
        self.__dict__.update(args)

class SoftTimer(object):
    deadline = None
    starttime = None
    checkcount = 0
    steps = None

    def __init__(self, timeout=None, **kwargs):
        self.set_timeout(timeout,from_now=True)
        self.steps=[]
        super(SoftTimer, self).__init__(**kwargs)


    def set_timeout(self, timeout, from_now=False):
        if from_now:
            self.starttime = time.time()
        if timeout and timeout > 0:
            self.deadline = self.starttime + timeout

    def check(self, step=None, logfn=None, **kwargs):
        self.checkcount +=1
        if step:
            self.steps.append(step)

        if logfn:
            logfn("Checking for timeout, %d time: %r", self.checkcount, step)

        t2 = time.time()
        if self.deadline and t2 > self.deadline:
            raise Timeout(tdiff=(t2 - self.starttime),
                          checkcount=self.checkcount,
                          steps = self.steps,
                          **kwargs)
