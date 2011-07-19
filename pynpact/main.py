import logging, os.path, tempfile, time, shutil
from optparse import OptionParser
from contextlib import contextmanager

from Bio import SeqIO

import util


#TODO: general next steps:

#start figuring out how to supply arguments to the various programs,
#and hash output filenames based on those args so that datestamp cache
#handles different invocations of a program.

#make paths less dependent on ~nathan

BINPATH=os.path.realpath(os.path.join(os.path.dirname(__file__),"../luciano-c/bin"))


def binfile(name) :
    return os.path.join(BINPATH, name)

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
    logger = logging.getLogger('pynpact.main')
    outputdir = None
    cleanup = True


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


    def run_extract(self,gene_descriptor="gene") :
        """Go through the genbank record pulling out gene names and locations
        $ extract MYCGE.gbk 0 gene 0 locus_tag > MYCGE.genes
"""

        outfilename,generate = self.derivative_filename(".extracted")

        if not generate :
            self.logger.debug("Skipped extraction")
            return outfilename

        with self.mkstemp_overwrite(outfilename) as outfile :
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
        return outfilename

    def run_CG(self) :
        """Do the CG ratio calculations.
$ CG MYCGE.gbk 1 580074 201 51 3 > MYCGE.CG200
"""
        outfilename,generate = self.derivative_filename("CG200")
        if generate :
            with self.mkstemp_overwrite(outfilename) as outfile :
                util.capturedCall([binfile("CG"), self.gbkfile, 1, len(self.seqrec), 201, 51, 3],
                                  stdout=outfile, logger=self.logger)
        return outfilename

    def run_atg(self) :
        # ../atg MYCGE.gbk > atg.txt
        pass

    def run_Allplots(self) :
        #need to generate AllPlots.def
        #/Users/luciano/src/Allplots 0 50000 5 1000 3 > XANCA.S-profiles.001.ps
        # fprintf(stderr,"\nUsage:\n\n  Allplots start interval lines [x-tics period_of_frames]\n");
        # fprintf(stderr,"\nstart                 Genome interval first base.");
        # fprintf(stderr,"\ninterval              Number of bases.\n");
        # fprintf(stderr,"\nlines                 Number of lines on page (one page).\n");
        # fprintf(stderr,"\nx-tics                Number of subdivisions.\n");
        # fprintf(stderr,"\nperiod_of_frame       Number of frames.\n");


        with self.mkdtemp() as dtemp :
            allplots_name = os.path.join(dtemp,"Allplots.def")
            self.logger.debug("writing %r", allplots_name)

            with open(allplots_name, 'w') as allplots :
                def ap_file(name) :
                    if name :
                        #absolute paths overran buffers in Allplots so
                        #use relative to keep names short.
                        path = os.path.relpath(name,dtemp)
                        allplots.write(path + "\n")
                    else :
                        allplots.write("None\n")

                #NB the "Plot Title" is disregarded, but that line should also contain the total number of bases
                allplots.write("%s %d\n" % ("FOOBAR", len(self.seqrec))) 	#Plot Title
                allplots.write("C+G\n")		#Nucleotide(s)_plotted (e.g.: C+G)
                allplots.write("Page 1\n")		#First-Page title
                allplots.write("page rest\n")	#Title of following pages

                ap_file(None)			#File_of_unbiased_CDSs
                ap_file(None)			#File_of_conserved_CDSs
                ap_file(None)			#File_of_new_CDSs
                ap_file(None)			#File_of_potential_new_CDSs
                ap_file(None)			#File_of_stretches_where_CG_is_asymmetric
                ap_file(self.run_extract())		#File_of_published_accepted_CDSs
                ap_file(None)			#File_of_published_rejected_CDSs
                ap_file(None)			#File_of_blocks_from_new_ORFs_as_cds
                ap_file(None)			#File_of_blocks_from_annotated_genes_as_cds
                ap_file(None)			#File_of_GeneMark_regions
                ap_file(None)			#File_of_G+C_coding_potential_regions
                ap_file(None)			#File_of_met_positions (e.g.:D 432)
                ap_file(None)			#File_of_stop_positions (e.g.:D 432)
                ap_file(None)			#File_of_tatabox_positions (e.g.:105.73 D 432 TATAAAAG)
                ap_file(None)			#File_of_capbox_positions
                ap_file(None)			#File_of_ccaatbox_positions
                ap_file(None)			#File_of_gcbox_positions
                ap_file(None)			#File_of_kozak_positions
                ap_file(None)			#File_of_palindrom_positions_and_size
                ap_file(self.run_CG())		#File_list_of_nucleotides_in_200bp windows.
                ap_file(None)			#File_list_of_nucleotides_in_100bp windows.

            i = 0
            ppage = 50000
            while i*ppage < len(self.seqrec) :
                outfilename,generate = self.derivative_filename(".%03d.ps" % (i+1,))
                if generate :
                    self.logger.debug("Starting Allplots for %r", os.path.basename(outfilename))
                    with util.mkstemp_overwrite(outfilename) as psout :
                        util.capturedCall([binfile("Allplots"), i*ppage, ppage, 5, 1000, 3], stdout=psout,
                                          logger=self.logger,cwd=dtemp)
                else :
                    self.logger.debug("Skipping Allplots for %r", os.path.basename(outfilename))
                i += 1


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
