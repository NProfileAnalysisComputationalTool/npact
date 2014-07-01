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
