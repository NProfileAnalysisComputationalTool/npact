import pytest
from path import path

from pynpact import parsing, genbank


def test_initial(gbkfile):
    result = parsing.initial(gbkfile)
    assert 'outputdir' in result
    assert result['filename'] == gbkfile
    assert 'basename' in result

def detect_format(config):
    parsing.detect_format(config)
    assert config['format']
    ddnafile = path(config['ddna'])
    assert ddnafile.exists()
    assert config['length'] == ddnafile.size
    assert ddnafile.open('r').read(14) == 'CAAATTGCGCTACA'

def test_detect_gbkformat(gbkconfig):
    detect_format(gbkconfig)

def test_detect_fnaformat(fnaconfig):
    detect_format(fnaconfig)

def test_detect_rawformat(rawconfig):
    detect_format(rawconfig)

def test_genbank(gbkfile):
    seqrec = genbank.parse_seq_rec(gbkfile)
    assert seqrec
    assert 0 < len(seqrec)


def test_isgbk(gbkconfig, fnaconfig):
    assert True == parsing.isgbk(gbkconfig)
    assert False == parsing.isgbk(fnaconfig)


def test_length_gbk(gbkconfig):
    assert 3204 == parsing.length(gbkconfig)
    assert 3204 == gbkconfig['length']


def test_length_fna(fnaconfig):
    assert 3204 == parsing.length(fnaconfig)
    assert 3204 == fnaconfig['length']


def test_derive_filename(gbkconfig):
    result = parsing.derive_filename(gbkconfig, 'asdf', 'foo')
    assert result.endswith('/NC_017123-asdf.foo')
