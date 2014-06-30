import logging
import os
import os.path
import hashlib
import tempfile
import shutil
import sys
from optparse import OptionParser
from contextlib import contextmanager
import math
from subprocess import PIPE

from __init__ import binfile, DATAPATH
import capproc
import prepare
import util

logger = logging.getLogger('pynpact.extract')
statuslog = logging.getLogger('pynpact.statuslog')


class Hasher(object):
    def __init__(self):
        self.state = hashlib.sha1()

    def hashdict(self, dict_):
        for k in sorted(dict_.keys()):
            val = dict_.get(k)
            if val is not None:
                self.state.update(k)
                self.state.update(str(val))

    def hashfiletime(self, filename):
        self.state.update(str(os.path.getmtime(filename)))

    def hash(self, str_):
        self.state.update(str_)

    def hexdigest(self):
        return self.state.hexdigest()


class BaseStep(object):
    def __init__(self, outputdir):
        self.outputdir = outputdir

    def derive_filename(self, basefile, hash, newext):
        "Build target filename based on identifying pieces"
        base = os.path.splitext(os.path.basename(basefile))[0]
        filename = '%s-%s.%s' % (base, hash, newext)
        return os.path.join(self.outputdir, filename)


class ExtractStep(BaseStep):
    keys = ['GeneDescriptorKey1', 'GeneDescriptorKey2',
            'GeneDescriptorSkip1', 'GeneDescriptorSkip2',
            'base_file']

    @staticmethod
    def get_hash(config):
        config = util.reducedict(config, ExtractStep.keys)
        h = Hasher()
        h.hashdict(config)
        h.hashfiletime(binfile('extract'))
        h.hashfiletime(config['base_file'])
        return config, h.hexdigest()

    @staticmethod
    def _extract(gbkfile, target_file, config):
        #TODO: check file existence here or in tqdaemon
        with open(target_file, 'w') as out:
            statuslog.debug(
                "Extracting genes in %s.", os.path.basename(gbkfile))
            cmd = [binfile("extract"), gbkfile,
                   config['GeneDescriptorSkip1'], config['GeneDescriptorKey1'],
                   config['GeneDescriptorSkip2'], config['GeneDescriptorKey2']]
            return capproc.capturedCall(
                cmd, check=True, logger=logger, stdout=out, stderr=sys.stderr)
        return target_file

    def run(self, config, executor):
        """Go through the genbank record pulling out gene names and locations
        $ extract MYCGE.gbk 0 gene 0 locus_tag > MYCGE.genes
        """
        config, hash = ExtractStep.get_hash(config)
        gbkfile = config['base_file']
        target_file = self.derive_filename(gbkfile, hash, 'genes')
        executor.queue(
            target_file, ExtractStep._extract, [gbkfile, target_file, config])
        return target_file


class AllplotsStep(BaseStep):
    keys = ['first_page_title', 'following_page_title',
            'length', 'start_base', 'end_base', 'period', 'bp_per_page',
            'nucleotides', 'alternate_colors', 'basename']
    file_keys = ['File_of_unbiased_CDSs',]  # and a whole bunch more

    @staticmethod
    def get_hash(config):
        h = Hasher()
        if config.get('run_extract'):
            h.hash(ExtractStep.get_hash(config)[1])
        keys = AllplotsStep.keys
        config = util.reducedict(config, keys)
        h.hashdict(config)
        h.hashfiletime(binfile('Allplots'))
        return config, h.hexdigest()

    def run(self, config, executor):
        h = Hasher()
        # If extract is to be run, kick it off while we have all the
        # configuration
        promises = []
        if config.get('run_extract'):
            filename = ExtractStep(self.outputdir).run(config, executor)
            config['File_of_published_accepted_CDSs'] = filename
            promises.append(filename)

        # Strip down to the config for this task only
        config = util.reducedict(config, AllplotsStep.keys + AllplotsStep.file_keys)


        bp_per_page = config['bp_per_page']
        start_base = config.pop('start_base')
        end_base = config.pop('end_base')
        h.hashdict(config)
        basehash = h.hexdigest()

        page_count = math.ceil(float(end_base - start_base) / bp_per_page)
        page_num = 1  # page number offset
        filenames = []
        # per-page loop
        while start_base < end_base:
            h = Hasher()
            h.hash(basehash)
            h.hash(page_num)
            h.hash(start_base)
            h.hash(end_base)
            config['start_base'] = start_base
            if start_base + bp_per_page < end_base:
                config['end_base'] = start_base + bp_per_page
            else:
                config['end_base'] = end_base
            psname = self.derive_filename(
                config['basename'], h.hexdigest(), 'ps')
            filenames.append(executor.queue(
                psname, AllplotsStep._ap,
                [psname, config, page_num, page_count]), after=promises)
            page_num += 1
            start_base += bp_per_page

        return filenames

    @staticmethod
    def _ap(target_file, pconfig, page_num, page_count):
        statuslog.info(
            "Generating %d pages of graphical output: %2d%%",
            page_count, round(100.0 * page_num / page_count))

        cmd = [binfile("Allplots"), "-q", "--stdin"]
        if pconfig.get('alternate_colors'):
            cmd.append("-C")

        # add the rest of the required args
        cmd += [pconfig['start_base'],
                pconfig['bp_per_page'],
                # TODO: move these into config
                5,     # lines on a page
                1000,  # Number of subdivisions
                pconfig['period'],
                pconfig['end_base']]
        with open(target_file, 'w') as out:
            with capproc.guardPopen(
                    cmd, stdin=PIPE, stdout=out, stderr=False,
                    logger=logger) as ap:
                # write the allplots.def file information through stdin.
                AllplotsStep.write_allplots_def(ap.stdin, pconfig, page_num)
                ap.stdin.close()
                ap.wait()
        return target_file

    @staticmethod
    def write_allplots_def(out, pconfig, page_num):
        def wl(line):
            "helper function for writing a line to the allplots file."
            if line:
                out.write(line)
            out.write('\n')

        # NB: the "Plot Title" is disregarded, but that line
        # should also contain the total number of bases
        wl("%s %d" % ("FOOBAR", pconfig['length']))

        # Nucleotide(s)_plotted (e.g.: C+G)
        wl('+'.join(pconfig['nucleotides']))
        # First-Page title
        wl(pconfig['first_page_title'].format(page_num))
        # Title of following pages
        wl(pconfig['following_page_title'].format(page_num))
        # write the rest of the filenames that allplots might read from
        for k in AllplotsStep.file_keys:
            f = pconfig.get(k)
            wl(f)
            wl(pconfig.get(f, "None"))
