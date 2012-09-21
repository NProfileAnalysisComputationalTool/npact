import logging, os.path
import re



import Bio.GenBank, Bio.GenBank.Scanner

#from Bio import SeqIO
#from Bio.GenBank import RecordParser,_RecordConsumer
#from Bio.GenBank.Scanner import GenBankScanner

import capproc
import util

logger = logging.getLogger(__name__)


class InvalidGBKException(Exception):
    """This class should only ever contain messages that are safe to present to the user."""
    pass


def reduce_genbank(gbkfile):
    """An attempt to create a version of the gbk file that has all the
    features but not the sequence in it, didn't end up being a
    siginificant savings.
"""
    def filterfun(outfile):
        with open(gbkfile,'r') as infile:
            for l in infile:
                outfile.write(l)
                if l.startswith("ORIGIN"):
                    outfile.write("//\n")
                    return

    return util.safe_produce_new(util.derivative_filename(gbkfile,".noseq", replace_ext=False),
                                 filterfun, logger=logger)



def open_parse_gb_rec(gbkfile, reduce_first=False):
    """Open the GenBank file using the underlying biopython libraries
    so we can get at the do_features keyword (False is generally quite
    a bit faster)

    Returns a the Bio.GenBank specialized record type.
"""
    if reduce_first:
        raise NotImplementedError("reduce_first option must be False for now")

    #rec = GenBankScanner().parse(open('NC_007912.gbk','r'), do_features=False) 
    #SeqIO.read(gbkfile,"genbank")

    with open(gbkfile,'r') as handle: 
        rp = Bio.GenBank.RecordParser()
        #rp._scanner = Bio.GenBank.Scanner.GenBankScanner()
        rp._consumer = Bio.GenBank._RecordConsumer()
        rp._scanner.feed(handle, rp._consumer, do_features=False)
        return rp._consumer.data
    
def open_parse_seq_rec(gbkfile, reduce_first=False, do_features=False):
    """Open the GenBank file using the underlying biopython libraries
    so we can get at the do_features keyword (False is generally quite
    a bit faster)

    Returns a SeqRecord object--the same as Bio.SeqIO.read(<file>,'genbank')
"""
    if reduce_first:
        raise NotImplementedError("reduce_first option must be False for now")

    logger.info("Parsing genbank file (features:%s): %r",
                do_features, gbkfile)

    #rec = GenBankScanner().parse(open('NC_007912.gbk','r'), do_features=False) 
    #SeqIO.read(gbkfile,"genbank")
    
    with open(gbkfile,'r') as handle: 
        rp = Bio.GenBank.FeatureParser()
        rp._consumer = Bio.GenBank._FeatureConsumer(rp.use_fuzziness,rp._cleaner)
        rp._scanner.feed(handle, rp._consumer, do_features=do_features)
        return rp._consumer.data

def make_seq_unknown(seq_record):
    seq_record.seq = Bio.Seq.UnknownSeq(len(seq_record), seq_record.seq.alphabet)



#TODO: This should probably be on disk so it can be shared amongst processes.
parse_cache = {}

def try_parse(abs_path, force=False):
    if not os.path.exists(abs_path):
        return None

    mtime = os.path.getmtime(abs_path)

    if not force:
        #cache_lines are None for miss, or (date, data) otherwise.
        cache_line = parse_cache.get(abs_path)
        if cache_line and cache_line[0] >= mtime:
            logger.debug("Cache hit for %s", abs_path)
            return cache_line[1]

    data = {'basename': os.path.basename(abs_path),
            'mtime': mtime,
            'filesize': util.pprint_bytes(os.path.getsize(abs_path)),
            }
    gbrec = None
    try:
        gbrec = open_parse_seq_rec(abs_path)
        if isinstance(gbrec.seq, Bio.Seq.UnknownSeq):
            raise InvalidGBKException("File contains no sequence data.")
        
        data['length'] = len(gbrec)
        data['id'] = gbrec.id
        data['date'] = gbrec.annotations.get('date')
        data['description'] = gbrec.description
        logger.debug("Pulled all data from the seqrec.")
    except InvalidGBKException:
        raise
    except:
        logger.debug("Failed parsing %r, trying regex search.", abs_path)
        
        match = re.search(r'(\d+) ?bp', open(abs_path).readline())
        if match :
            data['length'] = match.group(1)
        else :
            raise InvalidGBKException("Unable to find sequence length.")

    if 'length' in data:
        data['end_base'] = data['length']
        

    parse_cache[abs_path] = (mtime,data)
    return data

def default_config(abs_path):
    config = {
        'nucleotides':['c','g'],

        ##keys for extract.c
        'GeneDescriptorKey1': 'gene',
        'GeneDescriptorKey2': 'locus_tag',
        'GeneDescriptorSkip1': 0,
        'GeneDescriptorSkip2': 0,

        ##keys for CG:
        'window_size': 201,
        'step': 51,
        'period':3,

        #acgt_gamma:
        'skip_prediction': False,
        'significance': "0.001",

        ##allplots
        #http://docs.python.org/library/string.html#formatspec
        'first_page_title': None,
        'following_page_title': 'Page {0}',

        'bp_per_page': 50000,
        'start_base': 0,
        
        }


    if(abs_path):
        data = try_parse(abs_path)
        if data:
            config.update(data)
            config['first_page_title'] = config.get('description') or config.get('basename') or 'Page 1'
    return config


CONFIG_HELP_TEXT = {
    'start_base': "The base pair coordinate at which to start graphing.",
    'end_base': "The base pair coordinate at which to end graphing.",
    'length': "The length, in base pairs, of the genome being analyzed.",
    'first_page_title': "The title of the page containing the beginning of the genome.",
    'following_page_title': "The title of the pages after the first. Use {0} to get the page number.",
    'skip_prediction': "Should the acgt_gamma prediction be run? (click here to skip)",
    'significance': "What should the acgt_gamma prediction consider significant?",
    'nucleotides': "The bases to count the frequency of on the primary strand.",
    'alternate_colors': "Modifies graph colors from RGB to an alternative set that should be more easily distinguishable by people with certain color-blindness conditions.",
    }

significance_levels = ("0.01", "0.001", "0.0001")
significance_levels = zip(significance_levels, significance_levels)
