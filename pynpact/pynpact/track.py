import re


# see utils.js: ExtractParser too
class CDS (object):
    lineMatcher = re.compile(
        r'(?P<name>[^ ]+) (?P<complement>complement\()?'
        r'(?P<start_approx><)?(?P<start>\d+)\.\.'
        r'(?P<end_approx>>)?(?P<end>\d+)')

    def __init__(self, name=None, start=None, end=None,
                 complement=False, approximate=False, line=None,
                 start_approximate=False, end_approximate=False, **keys):
        (self.name, self.start, self.end, self.complement, self.approximate,
         self.line, self.start_approximate, self.end_approximate) = (
             name, start, end, complement, approximate, line,
             start_approximate, end_approximate)
        if line:
            self._loadline(line)

    def _loadline(self, line):
        m = CDS.lineMatcher.match(line)
        if(m):
            d = m.groupdict()
            self.name = d.get('name')
            self.start = int(d.get('start'))
            self.end = int(d.get('end'))
            self.complement = d.get('complement') is not None
            self.start_approximate = d.get('start_approx')
            self.end_approximate = d.get('end_approx')
            self.approximate = (d.get('start_approx') is not None
                                or d.get('end_approx') is not None)

    def __str__(self):
        if self.complement:
            return "{name} complement({start}..{end})\n" \
                .format(**self.__dict__)
        else:
            return "{name} {start}..{end}\n".format(**self.__dict__)


class Track(object):
    def __init__(self, name=None, data=[],
                 mycoplasma=False, filename=None, metadata={}, **keys):
        (self.name, self.data, self.mycoplasma, self.filename, self.metadata) \
            = (name, data, mycoplasma, filename, {})
        if len(self.data) > 0:
            i = 0
            for o in self.data:
                if not isinstance(o, CDS):
                    self.data[i] = CDS(**o)
                i += 1
        if filename and len(data) == 0:
            self._read()

    def _readcomment(self, line):
        (k, v) = line.split(':')
        k = k.strip(' \t\n\r#')
        v = v.strip(' \t\r\n')
        if k in self.__dict__.keys():
            self.__dict__[k] = v
        else:
            self.metadata[k] = v

    def _read(self):
        #print "Track._read %s" % len(self.data)
        self.data = []
        with open(self.filename) as f:
            for line in f.readlines():
                if(line.startswith('#')):
                    self._readcomment(line)
                else:
                    self.data.append(CDS(line=line))

    def write(self, filename=None):
        if not filename:
            filename = self.filename
        self.filename = filename
        with open(filename, 'w') as f:
            f.write("#name:%s\n" % self.name)
            f.write("#mycoplasma:%s\n" % self.mycoplasma)
            for (k, v) in self.metadata.items():
                f.write("#%s:%s\n" % (k, v))
            for o in self.data:
                f.write(str(o))

    def todict(self):
        d = self.__dict__.copy()
        for (i, o) in enumerate(self.data):
            d['data'][i] = o.__dict__.copy()
        return d

    def read(filename):
        return Track(filename=filename)
