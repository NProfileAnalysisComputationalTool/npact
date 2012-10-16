#!/usr/bin/env python
import logging, os.path, tempfile, shutil
from optparse import OptionParser
from contextlib import contextmanager
import math
from subprocess import PIPE

from __init__ import binfile, DATAPATH
import capproc
import prepare
import util
from softtimeout import SoftTimer, Timeout

logger = logging.getLogger(__name__)



#### One of the main ideas in here is that we name output products
#### with a hash of the configuration that went into them. That
#### filename then goes back in the master configuration, and might be
#### the input to a future computation (i.e. Allplots). The purpose of
#### this is to reuse the products when we can, but guarantee that
#### changes in the configuration cause the apropriate programs to
#### create new output.


class GenBankProcessor(object ):
    """Manipulate Genebank files in the pynpact project.

    Genbank files are documented at: http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html

    Takes GenBankProcessor(<gbkfile>, **kwargs)

    Where kwargs:

     * force: (default: False) recreate intermediate files even if
       they already exist.
     * logger: (default: None) what logger to use
     * outputdir: (default: directory of genebank file) What directory
       to create intermediate and PS files in.
     * cleanup: (default: False) Should the temporary files (not
       intermediate products) be deleted in case of errors
     * timeout: (default: None) Don't start a new processing step if
       this number of seconds has already elapsed

"""
    #
    gbkfile = None
    force = False
    logger = logger
    outputdir = None
    cleanup = True
    config = None
    timer = None

    def __init__(self, gbkfile = None, timeout=None, **kwargs):
        self.__dict__.update(kwargs)
        #for k in kwargs: setattr(self,k, kwargs[k])
        if not self.timer:
            self.timer = SoftTimer(timeout=timeout)

        self.parse(gbkfile)

    def parse(self, gbkfile):
        """Open up a GenBank file, trying to get the sequence record off of it."""

        self.gbkfile = gbkfile
        if not os.path.exists(gbkfile):
            raise Exception("Asked to parse nonexistant genebank file: %r", gbkfile)

        if self.config is None:
            self.config = prepare.default_config(self.gbkfile)

        if not self.outputdir:
            self.outputdir = os.path.dirname(os.path.realpath(self.gbkfile))



    ####################################################################
    ## Helper Functions

    def derivative_filename(self, part):
        """Build the filename of a derivative product of the gbk
        file. If the derivative file already exists return whether it
        is out of date"""

        if not part[0] == ".":
            part = "." + part

        base = os.path.splitext(os.path.basename(self.gbkfile))[0]
        outfilename = os.path.join(self.outputdir, base + part)

        return outfilename

    def mkstemp_overwrite(self, destination, **kwargs):
        "wrapper for util.mkstemp_overwrite; provides default dir, logger, cleanup args"
        kwargs.setdefault('dir', self.outputdir)
        kwargs.setdefault('logger', self.logger)
        kwargs.setdefault('cleanup', self.cleanup)
        return util.mkstemp_overwrite(destination, **kwargs)

    def safe_produce_new(self, filename, func, **kwargs):
        """A wrapper for util.safe_produce_new (create a new file in a
        multiprocess safe way) setting class defaults
        """
        kwargs.setdefault('dir', self.outputdir)
        kwargs.setdefault('logger', self.logger)
        kwargs.setdefault('cleanup', self.cleanup)
        kwargs.setdefault('force', self.force)
        return util.safe_produce_new(filename, func, **kwargs)

    @contextmanager
    def mkdtemp(self, ignore_errors=True, **kwargs):
        """A wrapper for tempfile.mkdtemp, sets defaults based on class variables"""
        kwargs.setdefault('dir', self.outputdir)
        path = tempfile.mkdtemp(**kwargs)
        try:
            yield path
        finally:
            if self.cleanup:
                self.logger.debug("Cleaning up mkdtemp %r", path)
                shutil.rmtree(path, ignore_errors=ignore_errors)

    #####

    def seqrec(self):
        if self._seqrec is None:
            self._seqrec = prepare.open_parse_seq_rec(self.gbkfile, do_features=True)
        return self._seqrec


    def run_extract(self):
        """Go through the genbank record pulling out gene names and locations
        $ extract MYCGE.gbk 0 gene 0 locus_tag > MYCGE.genes
        """
        config,hash = util.reducehashdict(self.config,
                                          ['GeneDescriptorKey1', 'GeneDescriptorKey2',
                                           'GeneDescriptorSkip1', 'GeneDescriptorSkip2',
                                           ])
        def _extract(out):
            self.timer.check("Extracting genes in %s." % os.path.basename(self.gbkfile))
            self.logger.info("starting extract for %s", self.gbkfile)
            cmd = [binfile("extract"), self.gbkfile,
                   config['GeneDescriptorSkip1'], config['GeneDescriptorKey1'],
                   config['GeneDescriptorSkip2'], config['GeneDescriptorKey2']]
            return capproc.capturedCall(cmd, stdout=out, logger=self.logger, check=True)
        filename = self.derivative_filename(".%s.genes" % hash)
        self.safe_produce_new(filename, _extract, dependencies=[self.gbkfile])
        self.config['File_of_published_accepted_CDSs'] = filename
        return filename

    def run_nprofile(self):
        """Do the CG ratio calculations.

        $ CG MYCGE.gbk 1 580074 201 51 3 > MYCGE.CG200
        """
        config,hash = util.reducehashdict(self.config,
                                          ['nucleotides', 'length','window_size','step','period'])
        outfilename = self.derivative_filename(".%s.nprofile" % hash)
        def _nprofile(out):
            self.timer.check("Calculating n-profile.")
            progargs = [binfile("nprofile"), '-b', ''.join(config["nucleotides"]),
                        self.gbkfile, 1, config['length'],
                        config['window_size'], config['step'], config['period']]
            return capproc.capturedCall(progargs, stdout=out, logger=self.logger, check=True)

        self.safe_produce_new(outfilename, _nprofile, dependencies=[self.gbkfile])
        self.config['File_list_of_nucleotides_in_200bp windows'] = outfilename
        return outfilename


    def acgt_gamma(self):
        "Run the acgt_gamma gene prediction program."
        if self.config.get('skip_prediction', False):
            self.logger.debug("Skipping prediction.")
            return

        assert os.path.exists(DATAPATH), \
               "Missing pynpact/data for acgt_gamma prediction. Expected at " + DATAPATH

        config,hash = util.reducehashdict(self.config,
                                          ['significance', 'GeneDescriptorSkip1'])

        #TODO: acgt actually takes the string of characters to skip, not the length.
        outdir = self.derivative_filename('.{0}.predict'.format(hash))
        if util.is_outofdate(outdir, self.gbkfile, DATAPATH):
            cmd = [binfile("acgt_gamma"), "-q", "-q", self.gbkfile]
            if config.has_key('significance'):
                cmd.append(config['significance'])


            with self.mkdtemp() as dtemp:
                self.timer.check("Predicting new gene locations.")
                self.logger.info("Starting prediction(sig:%s) program in %s",
                                 self.config.get('significance'),
                                 dtemp)
                capproc.capturedCall(cmd, cwd=dtemp, check=True,
                                     env={'BASE_DIR_THRESHOLD_TABLES':DATAPATH},
                                     logger=self.logger)
                #TODO, while deleting this is inconsistent (small
                #window): files will #be deleted out of outdir before
                #the directory is deleted and the rename comes a
                #moment later. During that time it will be
                #inconsistent.
                if os.path.exists(outdir):
                    self.logger.debug("Removing existing prediction output at %s", outdir)
                    shutil.rmtree(outdir)
                self.logger.debug("Renaming from %s to %s", dtemp, outdir)
                os.rename(dtemp, outdir)

        self.logger.debug("Adding prediction filenames to config dict.")
        #strip 4 characters off here b/c that's how acgt_gamma does it
        #at about lines 262-270
        j = lambda ext: os.path.join(outdir, os.path.basename(self.gbkfile)[:-4] + ext)
        self.config['File_of_new_CDSs'] = j(".newcds")
        self.config['File_of_published_rejected_CDSs'] = j(".modified")
        self.config['File_of_G+C_coding_potential_regions'] = j('.profiles')
        self.config['acgt_gamma_output'] = outdir




    ####################################################################
    ## Working with Allplots

    #the files in order that Allplots.def is expected to have.
    #NB: these are the keys into the config dictionary, the values are the filenames.
    AP_file_keys = ['File_of_unbiased_CDSs',
                    'File_of_conserved_CDSs',
                    'File_of_new_CDSs',
                    'File_of_published_rejected_CDSs',               #switched with "file_of_potential_new_CDs"
                    'File_of_stretches_where_CG_is_asymmetric',
                    'File_of_published_accepted_CDSs',
                    'File_of_potential_new_CDSs',                    #switched with "file_of_published_rejected_CDs"
                    'File_of_blocks_from_new_ORFs_as_cds',
                    'File_of_blocks_from_annotated_genes_as_cds',
                    'File_of_GeneMark_regions',
                    'File_of_G+C_coding_potential_regions',
                    'File_of_met_positions (e.g.:D 432)',
                    'File_of_stop_positions (e.g.:D 432)',
                    'File_of_tatabox_positions (e.g.:105.73 D 432 TATAAAAG)',
                    'File_of_capbox_positions',
                    'File_of_ccaatbox_positions',
                    'File_of_gcbox_positions',
                    'File_of_kozak_positions',
                    'File_of_palindrom_positions_and_size',
                    'File_list_of_nucleotides_in_200bp windows',
                    'File_list_of_nucleotides_in_100bp windows']

    def write_allplots_def(self, config, page_num, stream):
        "Writes out an Allplots.def for a single run through of allplots."


    def run_Allplots(self):
        """Actually invoke Allplots, once for each page we're generating.

Example invocation:
$ /Users/luciano/src/Allplots 0 50000 5 1000 3 > XANCA.S-profiles.001.ps

Documentation from Allplots.c
Usage:\n\n  Allplots start interval lines [x-tics period_of_frames]
start                 Genome interval first base.
interval              Number of bases.
lines                 Number of lines on page (one page).
x-tics                Number of subdivisions.
period_of_frame       Number of frames.

"""
        #figure out hashed filename of ps output.
        hashkeys = set(
            self.AP_file_keys + [
            'first_page_title', 'following_page_title',
            'length', 'start_base', 'end_base', 'period', 'bp_per_page',
            'nucleotides', 'alternate_colors']
            )

        #this hash is used down below for the final combined pages file
        config,hash = util.reducehashdict(self.config, hashkeys)
        #get the list of filenames this step depends on.
        dependencies = map(config.get, self.AP_file_keys)

        page_count = math.ceil((config['end_base'] - config['start_base']) * 1.0 / config['bp_per_page'])

        #build the individual ps page files.
        filenames = []

        #page number offset
        page_num = 1

        #the hash of individual pages shouldn't depend on the
        #end_base, so leave that out.  We keep updating the
        #start_base in *this* config (different than the general
        #one).
        pconfkeys = hashkeys.difference(set(["end_base"]))

        while (config['start_base'] < config['end_base']):
            if (config['start_base'] + config['bp_per_page']) >= config['end_base']:
                #we're on the last page: end base should be included in pconfkeys
                pconfkeys = hashkeys

            pconfig,phash = util.reducehashdict(config, pconfkeys)
            def _ap(psout):
                self.timer.check("Generating %d pages of graphical output: %2d%%" %
                                 (page_count, round(100.0 * page_num / page_count)))

                self.logger.debug("Starting Allplots page %d for %r",
                                  page_num, os.path.basename(self.gbkfile))

                cmd = [binfile("Allplots"), "-q", "--stdin",]
                if config.get('alternate_colors'):
                    cmd.append("-C")

                #add the rest of the required args
                cmd += [config['start_base'],
                        config['bp_per_page'],
                        #TODO: move these into config
                        5,    #lines on a page
                        1000, #Number of subdivisions
                        config['period'],
                        config['end_base']]

                with capproc.guardPopen(cmd, stdin=PIPE,
                                        stdout=psout,
                                        stderr=False,
                                        logger=self.logger) as ap:


                    self.logger.debug("Writing Allplots for %d", page_num)

                    #helper function for writing a line to the allplots file.
                    def ap_wl(str):
                        if str: ap.stdin.write(str)
                        ap.stdin.write('\n')

                    #NB the "Plot Title" is disregarded, but that line
                    #should also contain the total number of bases
                    ap_wl("%s %d" % ("FOOBAR", pconfig['length']))     #Plot Title
                    ap_wl('+'.join(pconfig['nucleotides']))            #Nucleotide(s)_plotted (e.g.: C+G)
                    ap_wl(pconfig['first_page_title'].format(page_num)) #First-Page title
                    ap_wl(pconfig['following_page_title'].format(page_num) ) #Title of following pages


                    for f in self.AP_file_keys:
                        ap_wl(pconfig.get(f,"None"))
                    ap.stdin.close()
                    ap.wait()

            psname = self.derivative_filename("%s.ps" % (phash))
            filenames.append(psname)
            self.safe_produce_new(psname, _ap, dependencies=dependencies)
            page_num += 1
            config['start_base'] += config['bp_per_page']



        def combine_ps_files(psout):
            #While combining, insert the special markers so that
            #it will appear correctly as many pages.
            self.timer.check("Combining pages into output.")
            self.logger.info("combining postscript files")
            first = True
            psout.write("%!PS-Adobe-2.0\n\n")
            psout.write("%%Pages: {0}\n\n".format(len(filenames)))
            idx = 1
            for psf in filenames:
                psout.write("%%Page: {0}\n".format(idx))
                with open(psf,'r') as infile:
                    infile.readline()
                    infile.readline()
                    psout.write(infile.readline())
                    for l in infile:
                        psout.write(l)
                idx += 1

        combined_ps_name = self.derivative_filename("%s.ps" %(hash,))
        self.config['combined_ps_name'] = combined_ps_name
        return self.safe_produce_new(combined_ps_name, combine_ps_files, dependencies=filenames)

    def run_ps2pdf(self):
        ps2pdf = util.which('ps2pdf')
        if ps2pdf:
            config,hash= util.reducehashdict(self.config, ['combined_ps_name'])
            pdf_filename = self.derivative_filename(".%s.pdf" % hash)
            def _ps2pdf(out):
                self.timer.check("Converting to PDF")
                self.logger.info("Converting to PDF")
                cmd = [ps2pdf, config['combined_ps_name'], '-']
                return capproc.capturedCall(cmd, stdout=out, logger=self.logger, check=True)
            self.safe_produce_new(pdf_filename, _ps2pdf, dependencies=[self.gbkfile])
            self.config['pdf_filename'] = pdf_filename
            return pdf_filename
        else:
            self.logger.info("Skipping ps2pdf translation, program not found.")
            return self.config['combined_ps_name']


    RUN_FNS=['run_extract','run_nprofile','acgt_gamma','run_Allplots','run_ps2pdf']

    def process(self):
        val = None
        for fn in self.RUN_FNS:
            val = getattr(self, fn)()
        return val


def process_all(path, config):
    gbp = GenBankProcessor(path, config=config)
    gbp.process()
    return gbp.config


if __name__ == '__main__':
    parser = OptionParser("""%prog <genebank file>""")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")
    parser.add_option("-f", "--force", action="store_true", dest="force",
                      help="Force recomputation of intermediate products")
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        exit(1)

    logging.basicConfig(level=(options.verbose and logging.DEBUG or logging.INFO),
                        format="%(asctime)s %(name)-10s %(levelname)-8s %(message)s",
                        datefmt='%H:%M:%S')

    gbkp = GenBankProcessor(args[0],force=options.force)
    gbkp.process()
