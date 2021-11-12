import tempfile
import pytest
import pyfastx
import hiscram.genome as hg
from hiscram.regions import Fragment


@pytest.fixture
def fasta():
    fname = tempfile.NamedTemporaryFile(delete=False).name
    with open(fname, "w") as fh:
        fh.write(">chr0\nAACCCGGGTT\n")
    yield pyfastx.Fasta(fname)


@pytest.fixture
def chrom():
    yield hg.Chromosome("chr0", 10)


@pytest.fixture
def genome():
    fname = tempfile.NamedTemporaryFile(delete=False).name
    with open(fname, "w") as fh:
        fh.write(">chr1\nAACCCAAACC\n")
        fh.write(">chr2\nGTGTGT\n")
    yield hg.Genome(pyfastx.Fasta(fname))


def test_chrom_len(chrom):
    assert len(chrom) == 10


def test_chrom_neg_coord():
    with pytest.raises(ValueError):
        hg.Chromosome("chr0", -10)


def test_chrom_boundaries(chrom):
    assert list(chrom.boundaries) == [0, 10]


def test_chrom_clean_frags(chrom):
    chrom.frags.append(Fragment("chr0", 10, 10))
    assert len(chrom.frags) == 2
    chrom.clean_frags()
    assert len(chrom.frags) == 1


def test_chrom_insert_middle(chrom):
    chrom.insert(2, Fragment("chr0", 3, 5))
    assert len(chrom.frags) == 3
    assert chrom.boundaries[1] == 2
    assert chrom.boundaries[2] == 4  # 10 + 11


def test_chrom_insert_left(chrom):
    chrom.insert(0, Fragment("chr0", 5, 9))
    assert len(chrom.frags) == 2
    assert list(chrom.boundaries) == [0, 4, 14]


def test_chrom_insert_right(chrom):
    chrom.insert(10, Fragment("chr0", 5, 9))
    assert len(chrom.frags) == 2
    assert list(chrom.boundaries) == [0, 10, 14]


def test_chrom_invert(chrom):
    chrom.invert(5, 6)
    assert len(chrom.frags) == 3
    assert [fr.is_reverse for fr in chrom.frags] == [False, True, False]
    assert list(chrom.boundaries) == [0, 5, 6, 10]


def test_chrom_invert_left(chrom):
    chrom.invert(0, 3)
    assert len(chrom.frags) == 2
    assert [fr.is_reverse for fr in chrom.frags] == [True, False]
    assert list(chrom.boundaries) == [0, 3, 10]


def test_chrom_delete(chrom):
    chrom.delete(2, 8)
    assert list(chrom.boundaries) == [0, 2, 4]


def test_chrom_get_seq(chrom, fasta):
    whole_chrom = [frag for frag in chrom.get_seq(fasta)]
    assert whole_chrom[0] == "AACCCGGGTT"


def test_chrom_get_seq_inv(chrom, fasta):
    chrom.invert(0, 3)
    frags_seq = [frag for frag in chrom.get_seq(fasta)]
    assert frags_seq == ["GTT", "CCGGGTT"]


def test_genome(genome):
    assert list(genome.chroms.keys()) == ["chr1", "chr2"]
    assert len(genome.chroms["chr1"]) == 10
    assert len(genome.chroms["chr2"]) == 6


def test_genome_insert(genome):
    genome.insert("chr1", 2, Fragment("chr2", 1, 4))
    assert len(genome.chroms["chr1"].frags) == 3


def test_genome_translocate(genome):
    genome.translocate("chr1", 3, Fragment("chr2", 0, 4))
    genome.get_seq


def test_genome_get_seq(genome):
    seqs = genome.get_seq()
    for chrom in genome.chroms.keys():
        assert len("".join(list(seqs[chrom]))) == len(genome.chroms[chrom])


def test_genome_get_breakpoints(genome):
    bps = genome.get_breakpoints()
    assert len(bps) == 0
