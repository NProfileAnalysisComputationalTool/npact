import json
import logging
import Bio
import Bio.GenBank
from Bio.GenBank.Scanner import GenBankScanner
from Bio.SeqFeature import ExactPosition


logger = logging.getLogger(__name__)


def parse_seq_rec(gbkfile, do_features=False):
    "Parse the gbkfile into a Bio.SeqRecord.SeqRecord"
    with open(gbkfile, 'r') as handle:
        # this is pretty much the exact code of Bio.Seq.read, except
        # we can access the do_features flag. Disabling feature
        # parsing is substantially faster.
        return GenBankScanner(debug=0).parse(handle, do_features=False)


# def parse_header(gbkfile):
#     """Parse just the header of the Genbank file

#     This is using internal methods from BioPython to do a relatively
#     minimal amount of work.

#     This doesn't work because the biggest thing we want out of it is
#     the length; but that's not there

#     """
#     from Bio.GenBank import _FeatureConsumer
#     from Bio.GenBank.utils import FeatureValueCleaner
#     from Bio.GenBank.Scanner import GenBankScanner
#     scanner = GenBankScanner()
#     # from Bio.GenBank.Scanner.InsdcScanner.parse
#     consumer = _FeatureConsumer(use_fuzziness=1,
#                                 feature_cleaner=FeatureValueCleaner())
#     with open(gbkfile) as f:
#         #from Bio.GenBank.Scanner.InsdcScanner.feed
#         scanner.set_handle(f)
#         if not scanner.find_start():
#             #Could not find (another) record
#             consumer.data = None
#             return False
#         scanner._feed_first_line(consumer, scanner.line)
#         scanner._feed_header_lines(consumer, scanner.parse_header())
#     return consumer.data




# def reduce_genbank(gbkfile):
#     """An attempt to create a version of the gbk file that has all the
#     features but not the sequence in it, didn't end up being a
#     siginificant savings.
# """
#     def filterfun(outfile):
#         with open(gbkfile, 'r') as infile:
#             for l in infile:
#                 outfile.write(l)
#                 if l.startswith("ORIGIN"):
#                     outfile.write("//\n")
#                     return

#     return util.safe_produce_new(
#         util.derivative_filename(gbkfile, ".noseq", replace_ext=False),
#         filterfun, logger=logger)


# def open_parse_gb_rec(gbkfile, reduce_first=False):
#     """Open the GenBank file using the underlying biopython libraries
#     so we can get at the do_features keyword (False is generally quite
#     a bit faster)

#     Returns a the Bio.GenBank specialized record type.
# """
#     if reduce_first:
#         raise NotImplementedError("reduce_first option must be False for now")

#     #rec = GenBankScanner().parse(open('NC_007912.gbk','r'), do_features=False)
#     #SeqIO.read(gbkfile,"genbank")

#     with open(gbkfile, 'r') as handle:
#         rp = Bio.GenBank.RecordParser()
#         #rp._scanner = Bio.GenBank.Scanner.GenBankScanner()
#         rp._consumer = Bio.GenBank._RecordConsumer()
#         rp._scanner.feed(handle, rp._consumer, do_features=False)
#         return rp._consumer.data


def open_parse_seq_rec(gbkfile, reduce_first=False, do_features=False):
    """Open the GenBank file using the underlying biopython libraries
    so we can get at the do_features keyword (False is generally quite
    a bit faster)

    Returns a SeqRecord object--the same as Bio.SeqIO.read(<file>,'genbank')"""
    if reduce_first:
        raise NotImplementedError("reduce_first option must be False for now")

    logger.info("Parsing genbank file (features:%s): %r",
                do_features, gbkfile)

    #rec = GenBankScanner().parse(open('NC_007912.gbk','r'), do_features=False)
    #SeqIO.read(gbkfile,"genbank")

    with open(gbkfile, 'r') as handle:
        rp = Bio.GenBank.FeatureParser()
        rp._consumer = Bio.GenBank._FeatureConsumer(
            rp.use_fuzziness, rp._cleaner)
        rp._scanner.feed(handle, rp._consumer, do_features=do_features)
        return rp._consumer.data


def _write_gbk_track_json(fh, dicts):
    fh.write('{"name":"Input CDSs", "type":"extracts", "active":true,\n')
    fh.write('"data":[\n')
    for r in dicts[:-1]:
        json.dump(r, fh)
        fh.write(',\n')
    json.dump(dicts[-1], fh)
    fh.write(']}')


def gbk_to_track_json(gbkfile, outfilename):
    rec = open_parse_seq_rec(gbkfile, do_features=True)
    rtn = []
    for feat in rec.features:
        if feat.type != 'CDS':
            continue
        d = feat.qualifiers.copy()
        d['start'] = feat.location.start.real
        d['end'] = feat.location.end.real
        d['approximate'] = not(isinstance(feat.location.start, ExactPosition)
                               and isinstance(feat.location.end, ExactPosition))
        d['complement'] = feat.location.strand == -1
        d['type'] = 'CDS'
        if not d.get('name'):
            d['name'] = d.get('locus_tag')
        for (k, v) in d.items():
            if isinstance(v, list) and len(v) == 1:
                d[k] = v[0]
        rtn.append(d)

    if isinstance(outfilename, file):
        _write_gbk_track_json(outfilename, rtn)
    else:
        with open(outfilename, 'w') as fh:
            _write_gbk_track_json(fh, rtn)
