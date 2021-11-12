# Structures and methods to represent SV breakpoints and fragments.
from __future__ import annotations
import copy
from typing import Optional, Tuple
from dataclasses import dataclass, field


@dataclass(order=True)
class Position:
    """A single position in the genome, defined by a
    chromosome, genomic position and optionally a sign."""

    chrom: str
    coord: int
    # True means positive strand
    sign: Optional[bool] = field(default=None, compare=False)

    def has_sign(self) -> bool:
        """Whether sign information is available."""
        return not self.sign is None


class Fragment:
    """A region representing a DNA sequence. Coordinates are 0-based and left-open.

    Attributes
    ----------

    chrom:
        Chromosome on which the fragment is located.
    start:
        Coordinate where the fragment starts
        on the chromosome (smallest).
    end: Coordinate where the fragment ends (largest).
    """

    def __init__(
        self, chrom: str, start: int, end: int, is_reverse: bool = False,
    ):
        if end < start:
            raise ValueError("end cannot be smaller than start.")
        if start < 0:
            raise ValueError("Coordinates cannot be negative.")
        self.chrom = chrom
        self.start = start
        self.end = end
        self.is_reverse = is_reverse

    def __len__(self) -> int:
        return self.end - self.start

    def middle(self) -> int:
        return (self.start + self.end) // 2

    def intersect(self, other: Fragment) -> int:
        """Return the length of intersection with another fragment.
        If there is no intersection, returns 0."""
        same_chrom = self.chrom == other.chrom
        separate = (self.start > other.end) or (self.end < other.start)
        if same_chrom and not separate:
            overlap_start = max(self.start, other.start)
            overlap_end = min(self.end, other.end)
            return overlap_end - overlap_start
        return 0

    def split(self, rel_pos: int) -> Tuple[Fragment, Fragment]:
        """Split the fragment at a positiion relative from the start."""
        left_frag, right_frag = copy.copy(self), copy.copy(self)
        left_frag.end = self.start + rel_pos
        right_frag.start = self.start + rel_pos
        return (left_frag, right_frag)

    def flip(self):
        """Change fragment's sign."""
        self.is_reverse = not self.is_reverse


class BreakPoint:
    """Defines a breakpoint in the genome.
    A breakpoint associates 2 different genomic positions
    from independent fragments.

    Attributes
    ----------

    pos1, pos2:
        The genomic positions connected by the breakpoint.
    frag1, frag2:
        Segments of DNA connected by the breakpoint
    """

    def __init__(
        self, pos1: Position, pos2: Position,
    ):

        # Properties are use to set/get fragments instead
        self._frag1 = None
        self._frag2 = None
        self.pos1, self.pos2 = pos1, pos2

    @property
    def signs(self) -> Optional[Tuple[str, str]]:
        """Return signs if available."""
        if self.has_signs():
            return self.pos1.sign, self.pos2.sign

    @staticmethod
    def _check_pos_at_frag_end(pos: Position, frag: Fragment) -> bool:
        """Check whether a position is at the correct end
        of a fragment."""
        right_chrom = pos.chrom == frag.chrom
        if pos.has_sign():
            if pos.sign:
                right_pos = pos.coord == frag.end
            else:
                right_pos = pos.coord == frag.start
        else:
            right_pos = pos.coord in [frag.start, frag.end]
        return right_chrom & right_pos

    @property
    def frag1(self):
        """Get the first fragment of the breakpoint"""
        return self._frag1

    @frag1.setter
    def frag1(self, frag: Fragment):
        if self._check_pos_at_frag_end(self.pos1, frag):
            self._frag1 = frag
        else:
            raise ValueError(
                "Attempted to set frag1 which does not match breakpoint.pos1."
            )

    @property
    def frag2(self, frag: Fragment):
        return self._frag2

    @frag2.setter
    def frag2(self, frag: Fragment):
        if self._check_pos_at_frag_end(self.pos1, frag):
            self._frag1 = frag
        else:
            raise ValueError(
                "Attempted to set frag2 which does not match breakpoint.pos2."
            )

    def has_signs(self) -> bool:
        """Whether sign information is available"""
        return self.pos1.has_sign() & self.pos2.has_sign()

    def has_frags(self) -> bool:
        """Whether fragments information is available."""
        return self.frag1 & self.frag2

    def connect(self, other) -> bool:
        """Whether two breakpoints could be connected."""
        if self.has_signs() and other.has_signs():
            compat_signs = self.pos2.sign != other.pos1.sign
        else:
            compat_signs = True
        compat_chroms = self.pos2.chrom == other.pos1.chrom
        return compat_chroms & compat_signs

    def overlap(self, other) -> int:
        """Amount of overlap between two breakpoint's fragments"""
        if self.pos2.sign != other.pos1.sign:
            return self.frag2.overlap(other.frag1)
        return 0

    def flip(self):
        """Ensure pos1 and frag1 correspond to the smallest position."""
        if self.pos1 > self.pos2:
            self.pos1, self.pos2 = self.pos2, self.pos1
