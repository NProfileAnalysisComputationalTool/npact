import os.path

from functools import wraps


class BaseStep(object):
    def __init__(self, config, outputdir, executor):
        self.config = config
        self.outputdir = outputdir
        self.executor = executor

    def derive_filename(self, hash, newext):
        "Build target filename based on identifying pieces"
        if newext[0] == '.':
            newext = newext[1:]
        base = os.path.splitext(self.config['basename'])[0]
        filename = '%s-%s.%s' % (base, hash, newext)
        return os.path.join(self.outputdir, filename)

    def replace_ext(self, base, newext):
        if newext[0] == '.':
            newext = newext[1:]
        return os.path.splitext(base)[0] + '.' + newext

    @classmethod
    def fromstep(klass, thisstep):
        return klass(thisstep.config, thisstep.outputdir, thisstep.executor)


class Task(object):
    """Small object to hold state for a function call to happen later.

    The point of this is to be a pickable closure looking thing.

    E.g.

        def adder(a,b):
            return a + b

        Task(adder, 1, 2)() == 3
    """
    func = None
    args = None
    kwargs = None

    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __call__(self):
        return self.func(*self.args, **self.kwargs)


def delay(fn):
    @wraps(fn)
    def wrapper(*args, **kwargs):
        return Task(fn, *args, **kwargs)
    return wrapper
