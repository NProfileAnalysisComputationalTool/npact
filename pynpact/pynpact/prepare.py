import logging, os.path, tempfile, time, shutil
from optparse import OptionParser
from contextlib import contextmanager

import Bio.GenBank, Bio.GenBank.Scanner

#from Bio import SeqIO
#from Bio.GenBank import RecordParser,_RecordConsumer
#from Bio.GenBank.Scanner import GenBankScanner

import util

logger = logging.getLogger(__name__)

def reduce_genbank(gbkfile) :
    def filterfun(outfile) :
        with open(gbkfile,'r') as infile:
            for l in infile :
                outfile.write(l)
                if l.startswith("ORIGIN") :
                    outfile.write("//\n")
                    return

    return util.safe_produce_new(gbkfile,".noseq", filterfun,
                                 replace_ext=False, logger=logger)



def open_parse_gb_rec(gbkfile, reduce_first=False) :
    if reduce_first :
        raise NotImplementedError("reduce_first option must be False for now")

    #rec = GenBankScanner().parse(open('NC_007912.gbk','r'), do_features=False) 
    #SeqIO.read(gbkfile,"genbank")

    with open(gbkfile,'r') as handle: 
        rp =Bio.GenBank.RecordParser()
        #rp._scanner = Bio.GenBank.Scanner.GenBankScanner()
        rp._consumer = Bio.GenBank._RecordConsumer()
        rp._scanner.feed(handle, rp._consumer, do_features=False)
        return rp._consumer.data
    
def open_parse_seq_rec(gbkfile, reduce_first=False) :
    if reduce_first :
        raise NotImplementedError("reduce_first option must be False for now")

    #rec = GenBankScanner().parse(open('NC_007912.gbk','r'), do_features=False) 
    #SeqIO.read(gbkfile,"genbank")

    with open(gbkfile,'r') as handle: 
        rp =Bio.GenBank.FeatureParser()

        rp._consumer = Bio.GenBank._FeatureConsumer(rp.use_fuzziness,rp._cleaner)
        rp._scanner.feed(handle, rp._consumer, do_features=False)
        return rp._consumer.data

def make_seq_unknown(seq_record) :
    seq_record.seq = Bio.Seq.UnknownSeq(len(seq_record),seq_record.seq.alphabet)
