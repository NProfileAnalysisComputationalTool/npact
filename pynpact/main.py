import logging
from Bio import SeqIO


def run_genes(gbk,record) :
    #this can somewhat be done with the SeqIO.features list, but not quite the same.
    pass

def parse(gbk) :
    return SeqIO.parse("input_files/MYCGE.gbk","genbank").next()

if __name__ == '__main__' :
    pass
