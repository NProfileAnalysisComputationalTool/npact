import os.path

from functools import wraps


class BaseStep(object):
    def __init__(self, config, outputdir, executor):
        self.config = config
        self.outputdir = outputdir
        self.executor = executor

    def derive_filename(self, hash, newext):
        "Build target filename based on identifying pieces"
        return derive_filename(self.config, hash, newext)

    @classmethod
    def fromstep(klass, thisstep):
        return klass(thisstep.config, thisstep.outputdir, thisstep.executor)


def derive_filename(config, hash, newext):
    "Build target filename based on identifying pieces"
    if newext[0] == '.':
        newext = newext[1:]
    base = os.path.splitext(config['basename'])[0]
    filename = '%s-%s.%s' % (base, hash, newext)
    return os.path.join(config['outputdir'], filename)
