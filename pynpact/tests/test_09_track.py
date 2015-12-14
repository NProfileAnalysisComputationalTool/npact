import pytest
from path import Path
from pynpact import track

TESTTRACK = "testdata/NC_007760.genes"
TESTCDS = ["ADEH_RS01085 234539..235147",
           "ADEH_RS01090 complement(235152..236021)",
           "ADEH_RS01115 complement(<240102..240860)"]


@pytest.fixture()
def trackfile():
    f = Path(__file__).dirname().joinpath(TESTTRACK)
    assert f.exists()
    return str(f)


def test_readcds():
    o = track.CDS(line=TESTCDS[0])
    assert o.start == 234539
    assert o.end == 235147
    assert o.name == 'ADEH_RS01085'
    assert not o.complement
    assert not o.approximate
    o = track.CDS(line=TESTCDS[2])
    assert o.start == 240102
    assert o.end == 240860
    assert o.name == 'ADEH_RS01115'
    assert o.complement
    assert o.approximate


def test_readtrack(trackfile):
    tr = track.Track(filename=trackfile)
    assert tr.name == 'Input CDSs'
    assert tr.mycoplasma
    assert tr.metadata.get('jim') == 'bob'
    assert len(tr.data) == 6
