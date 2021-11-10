import pytest
import hiscram.genome as hg
import hiscram.breakpoint as bp


@pytest.fixture
def chrom():
    yield hg.Chromosome("chr0", 100)


def test_chrom_len(chrom):
    assert len(chrom) == 100


def test_chrom_neg_coord():
    with pytest.raises(ValueError):
        hg.Chromosome("chr0", -10)


def test_chrom_boundaries(chrom):
    assert list(chrom.boundaries) == [0, 100]


def test_chrom_clean_frags(chrom):
    chrom.frags.append(bp.Fragment("chr0", 100, 100))
    assert len(chrom.frags) == 2
    chrom.clean_frags()
    assert len(chrom.frags) == 1


def test_chrom_insert_middle(chrom):
    chrom.insert(10, bp.Fragment("chr0", 33, 44))
    assert len(chrom.frags) == 3
    assert chrom.boundaries[1] == 10
    assert chrom.boundaries[2] == 21  # 10 + 11


def test_chrom_insert_left(chrom):
    chrom.insert(0, bp.Fragment("chr0", 5, 9))
    assert len(chrom.frags) == 2
    assert list(chrom.boundaries) == [0, 4, 104]


def test_chrom_insert_right(chrom):
    chrom.insert(100, bp.Fragment("chr0", 5, 9))
    assert len(chrom.frags) == 2
    assert list(chrom.boundaries) == [0, 100, 104]


def test_chrom_invert(chrom):
    chrom.invert(50, 60)
    assert len(chrom.frags) == 3
    assert [fr.is_reverse for fr in chrom.frags] == [False, True, False]
    assert list(chrom.boundaries) == [0, 50, 60, 100]


def test_chrom_delete(chrom):
    chrom.delete(12, 18)
    assert list(chrom.boundaries) == [0, 12, 94]


def test_chrom_get_seq():
    ...
