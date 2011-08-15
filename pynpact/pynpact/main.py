#!/usr/bin/env python
import logging, os.path, tempfile, time, shutil
from optparse import OptionParser
from contextlib import contextmanager

from Bio import SeqIO

from pynpact import binfile
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
        self.gbkfile = gbkfile

        if os.path.exists(gbkfile) :
            self.logger.info("Parsing genbank file: %r", gbkfile)
            #TODO: gbk can have multiple records in which case this
            #will err (i wasn't able to find one that did though)
            self.seqrec = SeqIO.read(gbkfile,"genbank")
        else :
            raise Exception("Asked to parse nonexistant genebank file: %r", gbkfile)

        if not self.outputdir :
            self.outputdir = os.path.dirname(os.path.realpath(self.gbkfile))


    def derivative_filename(self, part) :
        """Build the filename of a derivative product of the gbk
        file. If the derivative file already exists return whether it
        is out of date"""

        if not part[0] == "." :
            part = "." + part

        base = os.path.splitext(os.path.basename(self.gbkfile))[0]
        outfilename = os.path.join(self.outputdir, base + part)

        outofdate = self.force or (not os.path.exists(outfilename)) or \
                    os.path.getmtime(outfilename) < os.path.getmtime(self.gbkfile)

        return (outfilename, outofdate)

    def mkstemp_overwrite(self, destination, **kwargs) :
        kwargs.setdefault('dir', self.outputdir)
        kwargs.setdefault('logger', self.logger)
        kwargs.setdefault('cleanup', self.cleanup)
        return util.mkstemp_overwrite(destination, **kwargs)

    def safe_produce_new(self, ext_part, func, **kwargs) :
        kwargs.setdefault('dir', self.outputdir)
        kwargs.setdefault('logger', self.logger)
        kwargs.setdefault('cleanup', self.cleanup)
        kwargs.setdefault('force', self.force)
        return util.safe_produce_new(self.gbkfile, ext_part,func, **kwargs)

    @contextmanager
    def mkdtemp(self,**kwargs) :
        kwargs.setdefault('dir', self.outputdir)
        path = tempfile.mkdtemp(**kwargs)
        try :
            yield path
        finally :
            if self.cleanup :
                self.logger.debug("Cleaning up mkdtemp %r", path)
                shutil.rmtree(path,ignore_errors=True)

    def biopy_extract(self,gene_descriptor="gene") :

        def func(outfile):
            def print_feature(desc, strand, start, end) :
                location = "%s..%s" % (start +1, end)
                if f.strand == -1 :
                    location = "complement(" + location + ")"
                print >>outfile, "%-11s %s" % (desc, location )

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
        return self.safe_produce_new(".extracted", func)

    def original_extract(self) :
        func = lambda f: util.capturedCall([binfile("extract"), self.gbkfile, 0, "gene", 0, "locus_tag"],
                                           stdout=f, logger=self.logger)
        return self.safe_produce_new("genes", func)

    def run_extract(self) :
        """Go through the genbank record pulling out gene names and locations
        $ extract MYCGE.gbk 0 gene 0 locus_tag > MYCGE.genes
"""
        return self.biopy_extract()

    def run_CG(self) :
        """Do the CG ratio calculations.
$ CG MYCGE.gbk 1 580074 201 51 3 > MYCGE.CG200
"""
        progargs = [binfile("CG"), self.gbkfile, 1, len(self.seqrec), 201, 51, 3]
        func = lambda outfile: util.capturedCall(progargs, stdout=outfile, logger=self.logger)
        return self.safe_produce_new("CG200",func)

    def run_atg(self) :
        # ../atg MYCGE.gbk > atg.txt
        pass


    def write_allplots_def(self,allplots_name,page_num) :
        self.logger.debug("Writing %r", allplots_name)

        with open(allplots_name, 'w') as allplots :
            #need to generate AllPlots.def
            #/Users/luciano/src/Allplots 0 50000 5 1000 3 > XANCA.S-profiles.001.ps
            # fprintf(stderr,"\nUsage:\n\n  Allplots start interval lines [x-tics period_of_frames]\n");
            # fprintf(stderr,"\nstart                 Genome interval first base.");
            # fprintf(stderr,"\ninterval              Number of bases.\n");
            # fprintf(stderr,"\nlines                 Number of lines on page (one page).\n");
            # fprintf(stderr,"\nx-tics                Number of subdivisions.\n");
            # fprintf(stderr,"\nperiod_of_frame       Number of frames.\n");

            def ap_wl(str) :
                if str : allplots.write(str)
                allplots.write('\n')

            def ap_file(name=None) :
                #We may need to do other logic to calculate path, so separate function.
                ap_wl(name)

            #NB the "Plot Title" is disregarded, but that line should also contain the total number of bases
            ap_wl("%s %d" % ("FOOBAR", len(self.seqrec))) 	#Plot Title
            ap_wl("C+G")	      #Nucleotide(s)_plotted (e.g.: C+G)
            ap_wl(self.config.get('first_page_title' or "Page 1"))    #First-Page title
            ap_wl(self.config.get('following_page_title') or ("Page " + str(page_num))) #Title of following pages

            ap_file()                   #File_of_unbiased_CDSs
            ap_file()                   #File_of_conserved_CDSs
            ap_file()                   #File_of_new_CDSs
            ap_file()                   #File_of_potential_new_CDSs
            ap_file()                   #File_of_stretches_where_CG_is_asymmetric
            ap_file(self.run_extract()) #File_of_published_accepted_CDSs
            ap_file()                   #File_of_published_rejected_CDSs
            ap_file()                   #File_of_blocks_from_new_ORFs_as_cds
            ap_file()                   #File_of_blocks_from_annotated_genes_as_cds
            ap_file()                   #File_of_GeneMark_regions
            ap_file()                   #File_of_G+C_coding_potential_regions
            ap_file()                   #File_of_met_positions (e.g.:D 432)
            ap_file()                   #File_of_stop_positions (e.g.:D 432)
            ap_file()                   #File_of_tatabox_positions (e.g.:105.73 D 432 TATAAAAG)
            ap_file()                   #File_of_capbox_positions
            ap_file()                   #File_of_ccaatbox_positions
            ap_file()                   #File_of_gcbox_positions
            ap_file()                   #File_of_kozak_positions
            ap_file()                   #File_of_palindrom_positions_and_size
            ap_file(self.run_CG())	#File_list_of_nucleotides_in_200bp windows.
            ap_file()                   #File_list_of_nucleotides_in_100bp windows.


    def run_Allplots(self) :

        config,hash= util.reducehashdict(self.config, ['first_page_title','following_page_title'])

        with self.mkdtemp() as dtemp :
            i = 0
            ppage = 50000
            filenames = []
            while i*ppage < len(self.seqrec) :
                outfilename,generate = self.derivative_filename("%s.%03d.ps" % (hash,i+1))
                filenames.append(outfilename)
                if generate :
                    self.write_allplots_def(os.path.join(dtemp,"Allplots.def"), i+1)

                    self.logger.debug("Starting Allplots for %r", os.path.basename(outfilename))
                    with util.mkstemp_overwrite(outfilename) as psout :
                        util.capturedCall([binfile("Allplots"), i*ppage, ppage, 5, 1000, 3], stdout=psout,
                                          logger=self.logger,cwd=dtemp)
                else :
                    self.logger.debug("Skipping Allplots for %r", os.path.basename(outfilename))
                i += 1
            return filenames[0]


if __name__ == '__main__' :
    parser = OptionParser("""main.py <genebank file>""")
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
