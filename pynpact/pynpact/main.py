#!/usr/bin/env python
import logging, os.path, tempfile, time, shutil, re
from optparse import OptionParser
from contextlib import contextmanager

from Bio import SeqIO

from __init__ import binfile,DATAPATH
import util

logger = logging.getLogger(__name__)

#TODO: general next steps:

#start figuring out how to supply arguments to the various programs,
#and hash output filenames based on those args so that datestamp cache
#handles different invocations of a program.

#make paths less dependent on ~nathan


class GenBankProcessor(object) :
    """Manipulate Genebank files in the pynpact project.

    Genbank files are documented at: http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html

    Takes GenBankProcessor(<gbkfile>, **kwargs)

    Where kwargs:

     * force: (default: False) recreate intermediate files even if
       they already exist.
     * logger: (default: None) what logger to use
     * outputdir: (default: directory of genebank file) What directory
       to create intermediate and PS files in.
     * cleanup: (default: True) Should the temporary files (not
       intermediate products) be deleted in case of errors

"""
    #
    gbkfile = None
    force = False
    logger = logger
    outputdir = None
    cleanup = True
    config={}


    def __init__(self, gbkfile = None,**kwargs) :
        for k in kwargs : setattr(self,k, kwargs[k])

        if gbkfile :
            self.parse(gbkfile)

    def parse(self,gbkfile) :
        """Open up a GenBank file, trying to get the sequence record off of it."""
        self.gbkfile = gbkfile

        if os.path.exists(gbkfile) :
            self.logger.info("Parsing genbank file: %r", gbkfile)
            #TODO: gbk can have multiple records in which case this
            #will err (i wasn't able to find one that did though)
            try :
                self.seqrec = SeqIO.read(gbkfile,"genbank")
                self.config['length'] = len(self.seqrec)
            except:
                match = re.search(r'(\d+) ?bp', open(gbkfile).readline())
                if match :
                    self.config['length'] = match.group(1)
                else :
                    raise Exception("Unable to find sequence length on gbkfile: %r", gbkfile)

        else :
            raise Exception("Asked to parse nonexistant genebank file: %r", gbkfile)

        if not self.outputdir :
            self.outputdir = os.path.dirname(os.path.realpath(self.gbkfile))



    ####################################################################
    ## Helper Functions

    def derivative_filename(self, part) :
        """Build the filename of a derivative product of the gbk
        file. If the derivative file already exists return whether it
        is out of date"""

        if not part[0] == "." :
            part = "." + part

        base = os.path.splitext(os.path.basename(self.gbkfile))[0]
        outfilename = os.path.join(self.outputdir, base + part)

        return outfilename

    def mkstemp_overwrite(self, destination, **kwargs) :
        kwargs.setdefault('dir', self.outputdir)
        kwargs.setdefault('logger', self.logger)
        kwargs.setdefault('cleanup', self.cleanup)
        return util.mkstemp_overwrite(destination, **kwargs)

    def safe_produce_new(self, filename, func, **kwargs) :
        """A wrapper for util.safe_produce_new (create a new file in a
multiprocess safe way) setting class defaults
"""
        kwargs.setdefault('dir', self.outputdir)
        kwargs.setdefault('logger', self.logger)
        kwargs.setdefault('cleanup', self.cleanup)
        kwargs.setdefault('force', self.force)
        return util.safe_produce_new(filename, func, **kwargs)

    @contextmanager
    def mkdtemp(self,**kwargs) :
        """A wrapper for tempfile.mkdtemp, sets defaults based on class variables"""
        kwargs.setdefault('dir', self.outputdir)
        path = tempfile.mkdtemp(**kwargs)
        try :
            yield path
        finally :
            if self.cleanup :
                self.logger.debug("Cleaning up mkdtemp %r", path)
                shutil.rmtree(path,ignore_errors=True)

    def biopy_extract(self,gene_descriptor="gene") :
        """An implementation using the biopy library of the 'extract' functionality"""
        def func(outfile):
            def print_feature(desc, strand, start, end) :
                location = "%s..%s" % (start +1, end)
                if f.strand == -1 :
                    location = "complement(" + location + ")"
                print >>outfile, "%-11s %s" % (desc[:11], location )

            for f in self.seqrec.features :
                if f.type == 'CDS' :
                    desc = f.qualifiers.get(gene_descriptor)
                    #TODO: original stopped at the first space in the descripton
                    if desc :
                        if f.location_operator=="join" :
                            ex=1
                            for sf in f.sub_features :
                                #TODO: original one reversed the order
                                #here if it was on the complement
                                #strand
                                print_feature(desc[0] + "_E" + str(ex), sf.strand,
                                              sf.location.nofuzzy_start, sf.location.nofuzzy_end)
                                ex += 1
                        else :
                            print_feature(desc[0],f.strand,f.location.nofuzzy_start, f.location.nofuzzy_end)
        return self.safe_produce_new(self.derivative_filename("extracted"), func,
                                     dependencies=[self.gbkfile])
    
    def original_extract(self) :
        """Go through the genbank record pulling out gene names and locations
        $ extract MYCGE.gbk 0 gene 0 locus_tag > MYCGE.genes
"""
        config,hash = util.reducehashdict(self.config,
                                          ['GeneDescriptorKey1', 'GeneDescriptorKey2',
                                           'GeneDescriptorSkip1', 'GeneDescriptorSkip2',
                                           ])
        def func(f) :
            logger.info("starting extract for %s", self.gbkfile)
            cmd = [binfile("extract"), self.gbkfile,
                   config['GeneDescriptorSkip1'], config['GeneDescriptorKey1'],
                   config['GeneDescriptorSkip2'], config['GeneDescriptorKey2']]
            return util.capturedCall(cmd, stdout=f, logger=self.logger, check=True)
        filename = self.derivative_filename(".%s.genes" % hash)
        return self.safe_produce_new(filename, func, dependencies=[self.gbkfile])

    def run_extract(self) :
        """Go through the genbank record pulling out gene names and locations
"""
        return self.original_extract()

    def run_CG(self) :
        """Do the CG ratio calculations.
$ CG MYCGE.gbk 1 580074 201 51 3 > MYCGE.CG200
"""
        config,hash = util.reducehashdict(self.config,['length','window_size','step','period'])

        outfilename = self.derivative_filename(".%s.CG" % hash)
        progargs = [binfile("CG"), self.gbkfile, 1, config['length'],
                    config['window_size'], config['step'], config['period']]
        func = lambda out: util.capturedCall(progargs, stdout=out, logger=self.logger, check=True)
        return self.safe_produce_new(outfilename, func, dependencies=[self.gbkfile])


    def acgt_gamma(self):
        "Run the acgt_gamma gene prediction program."
        config,hash = util.reducehashdict(self.config,
                                          ['Significance', 'GeneDescriptorSkip1'])
        #TODO: acgt actually takes the string of characters to skip, not the length.
        outdir = self.derivative_filename('.{0}.predict'.format(hash))
        if util.is_outofdate(outdir, self.gbkfile):
            cmd = [binfile("acgt_gamma"), self.gbkfile]
            if config.has_key('Significance'):
                cmd.append(config['Significance'])
                
            
            with self.mkdtemp() as dtemp:
                self.logger.info("Starting prediction program in %s", dtemp)
                util.capturedCall(cmd, cwd=dtemp, check=True,
                                  env={'BASE_DIR_THRESHOLD_TABLES':DATAPATH},
                                  logger=self.logger)
                if os.path.exists(outdir):
                    self.loger.debug("Removing existing prediction output at %s", outdir)
                    shutil.rmtree(outdir)
                self.logger.debug("Renaming from %s to %s", dtemp, outdir)
                os.rename(dtemp,outdir)

        self.logger.debug("Adding prediction filenames to config dict.")
        #strip 4 characters off here b/c that's how acgt_gamma does it
        #at about lines 262-270
        gbkbase = os.path.filename(self.gbkfile)[:-4]
        j = lambda ext: os.path.join(outdir,gbkbase + ext)
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
                    'File_list_of_nucleotides_in_200bp windows.',
                    'File_list_of_nucleotides_in_100bp windows.']

    def write_allplots_def(self,config,allplots_name,page_num) :
        "Writes out an Allplots.def for a single run through of allplots."
        self.logger.debug("Writing %r", allplots_name)

        with open(allplots_name, 'w') as allplots :

            #helper function for writing a line to the allplots file.
            def ap_wl(str) :
                if str : allplots.write(str)
                allplots.write('\n')

            #NB the "Plot Title" is disregarded, but that line should also contain the total number of bases
            ap_wl("%s %d" % ("FOOBAR", config['length']))     #Plot Title
            ap_wl("C+G")                                      #Nucleotide(s)_plotted (e.g.: C+G)
            ap_wl(config['first_page_title'].format(page_num)) #First-Page title
            ap_wl(config['following_page_title'].format(page_num) ) #Title of following pages


            for f in self.AP_file_keys :
                ap_wl(config.get(f,"None"))
        return allplots_name


    def run_Allplots(self) :
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

        #do the dependent work.
        self.config['File_of_published_accepted_CDSs'] = self.run_extract()
        self.config['File_list_of_nucleotides_in_200bp windows.'] = self.run_CG()

        #figure out hashed filename of ps output.
        hashkeys = [ 'first_page_title', 'following_page_title', 'length'] + self.AP_file_keys
        config,hash= util.reducehashdict(self.config, hashkeys)

        #build the individual ps page files.
        filenames = []
        with self.mkdtemp() as dtemp :
            i = 0  #page number offset
            ppage = 50000
            while i*ppage < config['length'] :
                def dopage(psout) :

                    self.write_allplots_def(config, os.path.join(dtemp,"Allplots.def"), i+1)

                    self.logger.debug("Starting Allplots page %d for %r",
                                      i+1, os.path.basename(self.gbkfile))
                    
                    cmd = [binfile("Allplots"), i*ppage, ppage, 5, 1000, 3]
                    util.capturedCall(cmd, stdout=psout, stderr=False,
                                      logger=self.logger, cwd=dtemp, check=True)
                psname = self.derivative_filename("%s.%03d.ps" % (hash,i+1))
                filenames.append(psname)
                self.safe_produce_new(psname, dopage,
                                      dependencies=map(config.get,self.AP_file_keys))
                i += 1


            def combine_ps_files(psout) :
                #While combining, insert the special markers so that
                #it will appear correctly as many pages.
                self.logger.info("combining postscript files")
                first=True
                psout.write("%!PS-Adobe-2.0\n\n")
                psout.write("%%Pages: {0}\n\n".format(len(filenames)))
                idx=1
                for psf in filenames :
                    psout.write("%%Page: {0}\n".format(idx))
                    with open(psf,'r') as infile :
                        infile.readline()
                        infile.readline()
                        psout.write(infile.readline())
                        for l in infile :
                            psout.write(l)
                    idx += 1
            combined_ps_name = self.derivative_filename("%s.ps" %(hash,))
            return self.safe_produce_new(combined_ps_name, combine_ps_files, dependencies=filenames)


if __name__ == '__main__' :
    parser = OptionParser("""%prog <genebank file>""")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="Show more verbose log messages.")
    parser.add_option("-f","--force",action="store_true",dest="force",
                      help="Force recomputation of intermediate products")
    (options,args) = parser.parse_args()

    if len(args) != 1 :
        parser.print_help()
        exit(1)

    logging.basicConfig(level=(options.verbose and logging.DEBUG or logging.INFO),
                        format="%(asctime)s %(name)-10s %(levelname)-8s %(message)s",
                        datefmt='%H:%M:%S')

    gbkp = GenBankProcessor(args[0],force=options.force)
    gbkp.run_Allplots()
