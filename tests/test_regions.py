from copy import copy
import pytest
from hiscram.regions import BreakPoint, Fragment, Position

chroms = ["chr1", "chr1", "chr2", "chr2"]
coords = [100, 300, 1000, 10000]
signs = [False, True, False, None]

positions = [
    Position(chrom, coord, sign) for chrom, coord, sign in zip(chroms, coords, signs)
]

fragments = [Fragment(pos.chrom, pos.coord, pos.coord + 100) for pos in positions]

breakpoints = [BreakPoint(p1, p2) for p1, p2 in zip(positions[1:], positions[:-1])]


@pytest.mark.parametrize("pos", positions)
def test_position_has_sign(pos):
    assert pos.has_sign() == (False if pos.sign is None else True)


@pytest.mark.parametrize("frag", fragments)
def test_fragment_len(frag):
    assert len(frag) == 100


@pytest.mark.parametrize("frag", fragments)
def test_fragment_middle(frag):
    assert frag.middle() == frag.start + 50


@pytest.mark.parametrize("frag", fragments)
def test_fragment_intersect(frag):
    inter = frag.intersect(frag)
    assert inter == 100
    frag1 = copy(frag)
    frag1.start += 50
    inter = frag.intersect(frag1)
    assert inter == 50
    frag1.chrom = "chrXXX"
    inter = frag.intersect(frag1)
    assert inter == 0


@pytest.mark.parametrize("frag", fragments)
def test_fragment_split(frag):
    frag_l, frag_r = frag.split(0)
    assert len(frag_l) == 0
    assert len(frag_r) == 100
    frag_l, frag_r = frag.split(13)
    assert len(frag_l) == 13
    assert len(frag_r) == 87
    frag_l, frag_r = frag.split(100)
    assert len(frag_l) == 100
    assert len(frag_r) == 0


def test_breakpoint_signs():
    bp1 = BreakPoint(Position("chr1", 100, True), Position("chr1", 10000, False))
    bp2 = BreakPoint(Position("chr1", 300), Position("chr2", 100))
    assert bp1.has_signs()
    assert not bp2.has_signs()


@pytest.mark.parametrize("brp", breakpoints)
def test_breakpoint_frag1(brp):
    """Check if fragment setter and getter methods work as expected"""
    # Attempting to set illegal fragments
    with pytest.raises(ValueError):
        # wrong chrom
        brp.frag1 = Fragment("chrxx", 9999, 99999)
    with pytest.raises(ValueError):
        # wrong positions
        brp.frag1 = Fragment(brp.pos1.chrom, brp.pos1.coord + 1, brp.pos1.coord + 1)
    # Connected on the side matching sign
    if brp.pos1.has_sign():
        if brp.pos1.sign == True:
            brp.frag1 = Fragment(brp.pos1.chrom, brp.pos1.coord - 30, brp.pos1.coord)
        elif brp.pos1.sign == False:
            brp.frag1 = Fragment(brp.pos1.chrom, brp.pos1.coord, brp.pos1.coord + 1)
        else:
            brp.frag1 = Fragment(brp.pos1.chrom, brp.pos1.coord - 30, brp.pos1.coord)
            brp.frag2 = Fragment(brp.pos1.chrom, brp.pos1.coord, brp.pos1.coord + 1)
        assert brp.has_signs()


def test_breakpoint_connect():
    bp1 = breakpoints[0]
    bp2 = breakpoints[1]
    bp3 = breakpoints[2]
    assert bp1.connect(bp2)
    assert bp1.connect(bp3)


@pytest.mark.parametrize("brp", breakpoints)
def test_breakpoint_overlap(brp):
    ...


@pytest.mark.parametrize("brp", breakpoints)
def test_breakpoint_flip(brp):
    ...
