#!/usr/bin/env python
import logging, os.path, tempfile, time, shutil, re
from optparse import OptionParser
from contextlib import contextmanager

from Bio import SeqIO

from __init__ import binfile,DATAPATH
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

"""
    #
    gbkfile = None
    force = False
    logger = logger
    outputdir = None
    cleanup = False
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
        logger.info
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
    def mkdtemp(self,**kwargs):
        """A wrapper for tempfile.mkdtemp, sets defaults based on class variables"""
        kwargs.setdefault('dir', self.outputdir)
        path = tempfile.mkdtemp(**kwargs)
        try:
            yield path
        finally:
            if self.cleanup:
                self.logger.debug("Cleaning up mkdtemp %r", path)
                shutil.rmtree(path,ignore_errors=True)

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
        def thunk(out):
            self.timer.check("Extracting genes in %s." % os.path.basename(self.gbkfile))
            self.logger.info("starting extract for %s", self.gbkfile)
            cmd = [binfile("extract"), self.gbkfile,
                   config['GeneDescriptorSkip1'], config['GeneDescriptorKey1'],
                   config['GeneDescriptorSkip2'], config['GeneDescriptorKey2']]
            return util.capturedCall(cmd, stdout=out, logger=self.logger, check=True)
        filename = self.derivative_filename(".%s.genes" % hash)
        self.safe_produce_new(filename, thunk, dependencies=[self.gbkfile])
        self.config['File_of_published_accepted_CDSs'] = filename
        return filename

    def run_nprofile(self):
        """Do the CG ratio calculations.

        $ CG MYCGE.gbk 1 580074 201 51 3 > MYCGE.CG200
        """
        config,hash = util.reducehashdict(self.config,['length','window_size','step','period'])
        outfilename = self.derivative_filename(".%s.nprofile" % hash)
        def thunk(out):
            self.timer.check("Calculating n-profile.")
            progargs = [binfile("nprofile"), self.gbkfile, 1, config['length'],
                        config['window_size'], config['step'], config['period']]
            return util.capturedCall(progargs, stdout=out, logger=self.logger, check=True)
        
        self.safe_produce_new(outfilename, thunk, dependencies=[self.gbkfile])
        self.config['File_list_of_nucleotides_in_200bp windows'] = outfilename
        return outfilename
    

    def acgt_gamma(self):
        "Run the acgt_gamma gene prediction program."
        if not self.config.get('run_prediction', True):
            self.logger.debug("Skipping prediction.")
            return
        self.logger.debug("Starting prediction: %s", 
                          self.config.get('significance'))
        
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
                self.logger.info("Starting prediction program in %s", dtemp)
                util.capturedCall(cmd, cwd=dtemp, check=True,
                                  env={'BASE_DIR_THRESHOLD_TABLES':DATAPATH},
                                  logger=self.logger)
                #TODO, while deleting this is inconsistent: files will
                #be deleted out of outdir before the directory is
                #deleted and the rename comes a moment later. During
                #that time it will be inconsistent.
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




    ####################################################################
    ## Working with Allplots

    #the files in order that Allplots.def is expected to have.
    #NB: these are the keys into the config dictionary, the values are the filenames.
    AP_file_keys = ['File_of_unbiased_CDSs',
                    'File_of_conserved_CDSs',
                    'File_of_new_CDSs',
                    'File_of_potential_new_CDSs',
                    'File_of_stretches_where_CG_is_asymmetric',
                    'File_of_published_accepted_CDSs',
                    'File_of_published_rejected_CDSs',
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

    def write_allplots_def(self,config,allplots_name,page_num):
        "Writes out an Allplots.def for a single run through of allplots."
        self.logger.debug("Writing Allplots for %d %r", page_num, allplots_name)

        with open(allplots_name, 'w') as allplots:

            #helper function for writing a line to the allplots file.
            def ap_wl(str):
                if str: allplots.write(str)
                allplots.write('\n')

            #NB the "Plot Title" is disregarded, but that line should also contain the total number of bases
            ap_wl("%s %d" % ("FOOBAR", config['length']))     #Plot Title
            ap_wl("C+G")                                      #Nucleotide(s)_plotted (e.g.: C+G)
            ap_wl(config['first_page_title'].format(page_num)) #First-Page title
            ap_wl(config['following_page_title'].format(page_num) ) #Title of following pages


            for f in self.AP_file_keys:
                ap_wl(config.get(f,"None"))
        return allplots_name


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
        hashkeys = self.AP_file_keys + \
          ['first_page_title', 
           'following_page_title', 
           'length', 'start_page',
           'end_page', 'period', 'bp_per_page']
        config,hash= util.reducehashdict(self.config, hashkeys)

        #build the individual ps page files.
        filenames = []
        with self.mkdtemp() as dtemp:
            #page number offset
            page_num = config.get('start_page', 1)
            ppage = config['bp_per_page']

            #the hash of individual pages shouldn't depend on the
            #'start_page' and 'end_page' configuration, so we leave
            #that out.
            pconfkeys = set(hashkeys).difference(set(["start_page", "end_page"]))
            pconfig,phash = util.reducehashdict(config, pconfkeys)

            page_start = lambda page_num: ppage * (page_num - 1)

            while (page_start(page_num) < config['length']
                   and page_num <= config.get('end_page', 1000)):
                def thunk(psout):
                    self.timer.check("Generating page %d" % page_num)
                    self.write_allplots_def(pconfig, os.path.join(dtemp,"Allplots.def"), page_num)

                    self.logger.debug("Starting Allplots page %d for %r",
                                      page_num, os.path.basename(self.gbkfile))

                    cmd = [binfile("Allplots"), page_start(page_num), ppage, 5, 1000, config['period']]
                    util.capturedCall(cmd, stdout=psout, stderr=False,
                                      logger=self.logger, cwd=dtemp, check=True)
                psname = self.derivative_filename("%s.%03d.ps" % (phash, page_num))
                filenames.append(psname)
                self.safe_produce_new(psname, thunk,
                                      dependencies=map(config.get, self.AP_file_keys))
                page_num += 1


            def combine_ps_files(psout):
                #While combining, insert the special markers so that
                #it will appear correctly as many pages.
                self.timer.check("Combining pages into output.")
                self.logger.info("combining postscript files")
                first=True
                psout.write("%!PS-Adobe-2.0\n\n")
                psout.write("%%Pages: {0}\n\n".format(len(filenames)))
                idx=1
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
            return self.safe_produce_new(combined_ps_name, combine_ps_files, dependencies=filenames)

    RUN_FNS=['run_extract','run_nprofile','acgt_gamma','run_Allplots']

    def process(self):
        val = None
        for fn in self.RUN_FNS:
            val = getattr(self, fn)()
        return val





if __name__ == '__main__':
    parser = OptionParser("""%prog <genebank file>""")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")
    parser.add_option("-f","--force",action="store_true",dest="force",
                      help="Force recomputation of intermediate products")
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        exit(1)

    logging.basicConfig(level=(options.verbose and logging.DEBUG or logging.INFO),
                        format="%(asctime)s %(name)-10s %(levelname)-8s %(message)s",
                        datefmt='%H:%M:%S')

    gbkp = GenBankProcessor(args[0],force=options.force)
    gbkp.run_Allplots()
