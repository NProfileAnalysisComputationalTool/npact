import logging, os.path
from Bio import SeqIO

import util

BINPATH="/home/ACCELERATION/nathan/projects/spat/luciano-c/bin"
def binfile(name) :
    return os.path.join(BINPATH, name)

class GenBankProcessor(object) :
    #Genbank files are documented at: http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
    gbkfile = None
    force = False
    logger = logging

    def __init__(self, gbkfile = None) :
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


    def derivative_filename(self, part) :
        """Build the filename of a derivative product of the gbk
        file. If the derivative file already exists return whether it is out of date"""
        if not part[0] == "." :
            part = "." + part

        outfilename = os.path.splitext(self.gbkfile)[0] + part

        outofdate = self.force or (not os.path.exists(outfilename)) or \
                    os.path.getmtime(outfilename) < os.path.getmtime(self.gbkfile)

        return (outfilename, outofdate)


    def run_extract(self,gene_descriptor="gene") :
        """Go through the genbank record pulling out gene names and locations
$ extract MYCGE.gbk 0 gene 0 locus_tag > MYCGE.genes
"""


        outfilename,generate = self.derivative_filename(".extracted")

        if not generate :
            self.logger.debug("Skipped extraction")
            return outfilename

        with open(outfilename,'w') as outfile :
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
        outfilename,generate = self.derivative_filename("cg200ratio")
        if generate :
            with open(outfilename,'w') as outfile :
                util.capturedCall([binfile("CG"), self.gbkfile, 1, len(self.seqrec), 201, 51, 3],
                                  stdout=outfile,logger=self.logger)
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

        dirname = os.path.dirname(self.gbkfile)
        #TODO: temporary location for this so that it is thread/process safe
        allplots_name = os.path.join(dirname,"Allplots.def")
        self.logger.debug("writing %r", allplots_name)
        with open(allplots_name, 'w') as allplots :
            def ap_file(name) :
                if name :
                    #calculate the rel path.
                    path = os.path.relpath(name,dirname)
                    allplots.write(path + "\n")
                else :
                    allplots.write("None\n")

            #NB the "Plot Title" is disregarded, but that line should also contain the total number of bases
            allplots.write("%s %d\n" % ("FOOBAR",len(self.seqrec))) 	#Plot Title
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

        outfilename,generate = self.derivative_filename(".001.ps")
        with open(outfilename,'w') as psout :
            util.capturedCall([binfile("Allplots"), 0, 50000, 5, 1000, 3], stdout=psout,
                              logger=self.logger,cwd=dirname)

if __name__ == '__main__' :
    verbose = True
    logging.basicConfig(level=(verbose and logging.DEBUG or logging.INFO),
                        format="%(asctime)s %(name)-10s %(levelname)-8s %(message)s",
                        datefmt='%H:%M:%S')
    gbkp = GenBankProcessor()
    gbkp.parse("/home/ACCELERATION/nathan/projects/spat/input_files/MYCGE.gbk")
