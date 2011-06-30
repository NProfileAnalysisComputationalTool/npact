import logging, os.path
from Bio import SeqIO



class GenBankProcessor(object) :
    #Genbank files are documented at: http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
    gbkfile = None


    def __init__(self, gbkfile = None) :
        if gbkfile :
            self.parse(gbkfile)

    def run_extract(self,gene_descriptor="gene") :
        """Go through the genbank record pulling out gene names and locations
$ extract MYCGE.gbk 0 gene 0 locus_tag > MYCGE.genes
"""

        #this can somewhat be done with the SeqIO.features list, but not quite the same.


        outfilename = os.path.splitext(self.gbkfile)[0] + ".extracted"

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


    def parse(self,gbkfile) :
        self.gbkfile = gbkfile
            #TODO: gbk can have multiple records in which case this will err (i wasn't able to find one that did though)
        self.seqrec = SeqIO.read(gbkfile,"genbank")


    def run_CG(self) :
        """Do the CG ratio calculations.
$ CG MYCGE.gbk 1 580074 201 51 3 > MYCGE.CG200
"""
        pass

    def run_atg(self) :
        # ../atg MYCGE.gbk > atg.txt
        pass

    def run_Alllots(self) :
        #need to generate AllPlots.def
        #/Users/luciano/src/Allplots 0 50000 5 1000 3 > XANCA.S-profiles.001.ps
        pass


if __name__ == '__main__' :
    gbkp = GenBankProcessor()
    gbkp.parse("/home/ACCELERATION/nathan/projects/spat/input_files/MYCGE.gbk")
