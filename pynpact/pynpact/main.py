#!/usr/bin/env python
import logging
import os, os.path
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

logger = logging.getLogger('pynpact')



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

"""
    #
    gbkfile     = None
    force       = False
    logger      = logger #logger is for internal debugging purposes
    statuslog   = None   #logger for printing external messages to end users 
    outputdir   = None
    cleanup     = True
    config      = None


    def __init__(self, gbkfile = None, **kwargs):
        self.__dict__.update(kwargs)

        if self.config and self.config.get('raiseerror'):
            raise Exception("You asked me to.")

        #we want this to just print raw messages to stderr; they
        #will be piped onto the webpage.
        self.statuslog = logging.getLogger('pynpact.statuslog')
        self.statuslog.setLevel(logging.DEBUG)
        
        sh = logging.StreamHandler()
        sh.setFormatter(logging.Formatter('%(message)s'))
        self.statuslog.addHandler(sh)

        if self.config and self.config.get('force'):
            self.statuslog.debug("Forcing recalculation")
            self.force = self.config.get('force')
        
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
            self.statuslog.debug("Extracting genes in %s.", os.path.basename(self.gbkfile))
            cmd = [binfile("extract"), self.gbkfile,
                   config['GeneDescriptorSkip1'], config['GeneDescriptorKey1'],
                   config['GeneDescriptorSkip2'], config['GeneDescriptorKey2']]
            return capproc.capturedCall(cmd, 
                                        stdout=out, stderr=sys.stderr,
                                        logger=self.logger, check=True)
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
            self.statuslog.info("Calculating n-profile.")
            progargs = [binfile("nprofile"), '-b', ''.join(config["nucleotides"]),
                        self.gbkfile, 1, config['length'],
                        config['window_size'], config['step'], config['period']]
            return capproc.capturedCall(progargs, 
                                        stdout=out, stderr=sys.stderr,
                                        logger=self.logger, check=True)

        self.safe_produce_new(outfilename, _nprofile, dependencies=[self.gbkfile])
        self.config['File_list_of_nucleotides_in_200bp windows'] = outfilename
        return outfilename


    def acgt_gamma(self):
        "Identifying ORFs with significant 3-base periodicities."
        if self.config.get('skip_prediction', False):
            self.statuslog("Skipping ORF identification.")
            return

        assert os.path.exists(DATAPATH), \
               "Missing pynpact/data for acgt_gamma prediction. Expected at " + DATAPATH

        config,hash = util.reducehashdict(self.config,
                                          ['significance', 'GeneDescriptorSkip1'])

        #TODO: acgt actually takes the string of characters to skip, not the length.
        outdir = self.derivative_filename('.{0}.predict'.format(hash))
        if self.force or util.is_outofdate(outdir, self.gbkfile, DATAPATH):
            cmd = [binfile("acgt_gamma"), "-q", self.gbkfile]
            if config.has_key('significance'):
                cmd.append(config['significance'])


            with self.mkdtemp() as dtemp:
                self.statuslog.info("Identifying ORFs with 3-base compositional periodicity with significance: %s.",
                                    self.config.get('significance'))
                self.logger.debug("Starting prediction program in %s", dtemp)
                capproc.capturedCall(cmd, cwd=dtemp, check=True,
                                     env={'BASE_DIR_THRESHOLD_TABLES': DATAPATH},
                                     stderr=sys.stderr,
                                     logger=self.logger)
                #TODO, while deleting this is inconsistent (small
                #window): files will be deleted out of outdir before
                #the directory is deleted and the rename comes a
                #moment later. During that time it will be
                #inconsistent.
                if os.path.exists(outdir):
                    self.logger.debug("Removing existing output at %s", outdir)
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

        page_count = math.ceil(float(config['end_base'] - config['start_base'])
                               / config['bp_per_page'])

        #build the individual ps page files.
        filenames = []

        #page number offset
        page_num = 1

        #the hash of individual pages shouldn't depend on the
        #end_base, so leave that out.  We keep updating the start_base
        #in *this* config (different than the general one).
        pconfkeys = hashkeys.difference(set(["end_base"]))

        while (config['start_base'] < config['end_base']):
            if (config['start_base'] + config['bp_per_page']) >= config['end_base']:
                #we're on the last page: include end_base for hash
                #generation as the end_base might abridge the page's
                #output
                pconfkeys = hashkeys

            pconfig,phash = util.reducehashdict(config, pconfkeys)
            def _ap(psout):
                self.statuslog.info("Generating %d pages of graphical output: %2d%%",
                                    page_count, round(100.0 * page_num / page_count))

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
                    #write the allplots.def file information through stdin.
                    def ap_wl(line):
                        "helper function for writing a line to the allplots file."
                        if line: ap.stdin.write(line)
                        ap.stdin.write('\n')

                    #NB the "Plot Title" is disregarded, but that line
                    #should also contain the total number of bases
                    ap_wl("%s %d" % ("FOOBAR", pconfig['length']))           #Plot Title
                    ap_wl('+'.join(pconfig['nucleotides']))                  #Nucleotide(s)_plotted (e.g.: C+G)
                    ap_wl(pconfig['first_page_title'].format(page_num))      #First-Page title
                    ap_wl(pconfig['following_page_title'].format(page_num) ) #Title of following pages
                    for f in self.AP_file_keys:
                        ap_wl(pconfig.get(f, "None"))
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
            self.statuslog.info("Combining pages into output.")
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
                self.statuslog.info("Converting to PDF")
                cmd = [ps2pdf, config['combined_ps_name'], '-']
                capproc.capturedCall(cmd, stdout=out, logger=self.logger, check=True)
                self.statuslog.info("Finished PDF conversion")

            self.safe_produce_new(pdf_filename, _ps2pdf, dependencies=[self.gbkfile])
            self.config['pdf_filename'] = pdf_filename
            return pdf_filename
        else:
            self.logger.warn("Skipping ps2pdf translation, program not found.")
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

    original = args[0]
    gbkfile = os.path.realpath(original)
    logging.basicConfig(level=(options.verbose and logging.DEBUG or logging.INFO),
                        format="%(asctime)s %(name)-10s %(levelname)-8s %(message)s",
                        datefmt='%H:%M:%S')

    gbkp = GenBankProcessor(gbkfile, force=options.force)
    outputfile = gbkp.process()
    logging.info("Finished processing %r", original)
    logging.info("See output at %r", outputfile)
