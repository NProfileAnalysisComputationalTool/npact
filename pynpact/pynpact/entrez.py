import os, os.path
import datetime, time
import logging

import Bio.Entrez
import util

Bio.Entrez.email = "nathan@acceleration.net"
Bio.Entrez.tool = "Biopython->pynpact"


MAX_TO_SUMMARIZE = 10

logger = logging.getLogger(__name__)

class TooManyResponsesException(Exception):
    term=None
    summaries=None
    count=None
    WebEnv=None
    QueryKey=None
    def __init__(self,term=None,count=None,Webenv=None,QueryKey=None):
        term=self.term
        summaries=self.summaries
        count=self.count
        WebEnv=self.WebEnv
        QueryKey=self.QueryKey

class EntrezSession(object):
    db = 'genome'
    WebEnv = None
    QueryKey = None
    lib_path = None
    summaries = None

    def __init__(self, lib_path, **kwargs):
        self.lib_path = lib_path
        self.__dict__.update(kwargs)

    def reset(self):
        self.QueryKey = None
        self.WebEnv = None
        self.result_count = None

    def has_session(self):
        return self.QueryKey and len(self.QueryKey) and \
               self.WebEnv and len(self.WebEnv)

    def search(self, term):
        logger.info("Starting Entrez query for %r, session=%s", term, self.has_session())
        resp = Bio.Entrez.read(Bio.Entrez.esearch(db=self.db, term=term,
                                                  usehistory=True,
                                                  query_key=self.QueryKey,
                                                  webenv=self.WebEnv))
        self.QueryKey = resp['QueryKey']
        self.WebEnv = resp['WebEnv']
        self.result_count = int(resp['Count'])
        logger.debug("Got back %s results.", self.result_count)
        return resp

    def summarize(self):
        if not self.summaries:
            logger.info("Summarizing from %s, %s", self.WebEnv, self.QueryKey)
            self.summaries = Bio.Entrez.read(Bio.Entrez.esummary(db=self.db,
                                                                 webenv=self.WebEnv,
                                                                 query_key=self.QueryKey))
        return self.summaries

    def fetch(self, summary=None, filename=None):
        if not summary:
            if self.result_count == 1:
                summary = self.summarize()[0]
            else :
                raise TooManyResponsesException()
        id = summary['Id']
        #this appears to be the RefSeq Id.
        caption = summary['Caption']
        logger.info("Starting fetch of Id: %s, RefSeq ID: %s", id, caption)

        if filename is None:
            filename =  os.path.join(self.lib_path, caption + ".gbk")

        update_date = datetime.datetime(*map(int,summary['UpdateDate'].split('/')))

        if (not os.path.exists(filename)) or \
               datetime.datetime.fromtimestamp(os.path.getmtime(filename)) < update_date:
            #file should be downloaded.
            net_handle = Bio.Entrez.efetch(db=self.db, rettype='gb', id=id)
            logger.debug("Streaming handle to file.")
            with util.mkstemp_overwrite(filename,logger=logger) as f:
                bytes = util.stream_to_file(net_handle,f)
                logger.info("Saved %s to %s.", util.pprint_bytes(bytes), filename)
        return filename

    def fetch_id(self,id):
        logger.info("Starting fetch of Id: %s", id)
        summaries = Bio.Entrez.read(Bio.Entrez.esummary(db=self.db,id=id))
        return self.fetch(summary=summaries[0])

    def to_url(self,term):
        """Convert the query we've done to a url that will load ncbi's site."""
        fmt = "http://www.ncbi.nlm.nih.gov/sites/entrez?db={0}&term={1}"
        return fmt.format(self.db,term)

    def to_session_url(self):
        """Convert the query we've done to a url that will load ncbi's site."""
        fmt = "http://www.ncbi.nlm.nih.gov/sites/entrez?db={0}&cmd=HistorySearch&querykey={1}&tab=&WebEnv={2}"
        return fmt.format(self.db,self.QueryKey,self.WebEnv)


ENTREZ_CACHE={}
class CachedEntrezSession(EntrezSession):
    def search(self, term):
        self.term = term
        self.summaries = ENTREZ_CACHE.get(term)
        if self.summaries:
            logger.debug('CachedEntrezSession hit for %s', term)
            self.result_count = len(self.summaries)
        else:
            result = super(CachedEntrezSession, self).search(term)

    def summarize(self):
        if not self.summaries:
            ENTREZ_CACHE[self.term] = super(CachedEntrezSession, self).summarize()
        return self.summaries
        
        
    

# [{'Status': 'Completed',
#   'Comment': '  ',
#   'Caption': 'NC_014248',
#   'Title': "'Nostoc azollae' 0708 chromosome, complete genome",
#   'CreateDate': '2010/06/16',
#   'Extra': 'gi|298489614|ref|NC_014248.1||gnl|NCBI_GENOMES|26219[26219]',
#   'TaxId': 551115,
#   'ReplacedBy': '',
#   u'Item': [],
#   'Length': 5354700,
#   'Flags': 256,
#   'UpdateDate': '2011/04/06',
#   u'Id': '26219',
#   'Gi': 2621}
#  ]
