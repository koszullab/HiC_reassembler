"""Structures and methods to represent SV breakpoints and fragments.
This submodule defines a BreakPoint class representing connections between positions
on separate fragments. Each breakpoint consists of two positions. Each of position is
located on a fragment. Positions have a sign representing what end of a fragment
they're on.
"""
#
from __future__ import annotations
import copy
from typing import Optional, Tuple
from dataclasses import dataclass, field


@dataclass(order=True)
class Position:
    """A single position in the genome, defined by a
    chromosome, genomic position and optionally a sign.

    Attributes
    ----------

    chrom:
        The chromosome name.
    coord:
        The 0-based coordinate on the chromosome.
    sign:
        Whether the position is on the 5' (-) or 3' (+) side.

    """

    chrom: str
    coord: int
    #  None means unknown, True means 3'
    sign: Optional[bool] = field(default=None, compare=False)

    def __repr__(self):
        sign_symbol = {None: "", True: "+", False: "-"}
        return f"{self.chrom}:{self.coord}:{sign_symbol[self.sign]}"

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
        self,
        chrom: str,
        start: int,
        end: int,
        is_reverse: bool = False,
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

    def __repr__(self) -> str:
        sign = "-" if self.is_reverse else "+"
        return f"{self.chrom}:{self.start}-{self.end}:{sign}"

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

    def merge(self, other: Fragment) -> Optional[Fragment]:
        """Merge two fragments with adjacent genomic positions. If fragments
        cannot be merged (e.g. because they are not adjacent), returns None."""
        compat_chroms = self.chrom == other.chrom
        compat_signs = self.is_reverse != other.is_reverse
        compat_coords = (self.start == other.end) or (other.start == self.end)
        if compat_chroms and compat_signs and compat_coords:
            start = min(self.start, other.start)
            end = min(self.end, other.end)
            frag = Fragment(self.chrom, start, end, self.is_reverse)
        else:
            frag = None
        return frag


class BreakPoint:
    """Defines a breakpoint in the genome.
    A breakpoint associates 2 different genomic positions from independent
    fragments. pos1 and frag1 represent the 5' side of the breakpoint (I think ?).

    Attributes
    ----------

    pos1, pos2:
        The genomic positions connected by the breakpoint.
    frag1, frag2:
        Segments of DNA connected by the breakpoint
    """

    def __init__(
        self,
        pos1: Position,
        pos2: Position,
    ):

        # Properties are use to set/get fragments instead
        self._frag1 = None
        self._frag2 = None
        self.pos1, self.pos2 = pos1, pos2

    def __repr__(self) -> str:
        if self.has_frags():
            p1 = self.frag1
            p2 = self.frag2
        else:
            p1 = self.pos1
            p2 = self.pos2
        return f"{str(p1)}|{str(p2)}"

    @property
    def signs(self) -> Optional[Tuple[str, str]]:
        """Return signs if available."""
        if self.has_signs():
            sign1 = "+" if self.pos1.sign else "-"
            sign2 = "+" if self.pos2.sign else "-"
            return sign1, sign2

    @staticmethod
    def _check_pos_at_frag_end(pos: Position, frag: Fragment) -> bool:
        """Check whether a position is at the correct end
        of a fragment."""
        right_chrom = pos.chrom == frag.chrom
        if pos.has_sign():
            if pos.sign:
                if frag.is_reverse:
                    right_pos = pos.coord == frag.start
                else:
                    right_pos = pos.coord == frag.end
            else:
                if frag.is_reverse:
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
                f"Attempted to set frag1 to {frag} which does not match "
                f"breakpoint.pos1 of {self.pos1}."
            )

    @property
    def frag2(self):
        return self._frag2

    @frag2.setter
    def frag2(self, frag: Fragment):
        if self._check_pos_at_frag_end(self.pos2, frag):
            self._frag2 = frag
        else:
            raise ValueError(
                f"Attempted to set frag2 to {frag}, which does not match "
                f"breakpoint.pos2 of {self.pos2}."
            )

    def has_signs(self) -> bool:
        """Whether sign information is available"""
        return self.pos1.has_sign() & self.pos2.has_sign()

    def has_frags(self) -> bool:
        """Whether fragments information is available."""
        return (self.frag1 is not None) & (self.frag2 is not None)

    def can_connect(
        self, other: BreakPoint, min_intersect: int = 1000
    ) -> bool:
        """Whether two breakpoints could be connected (i.e. whether they could
        share a fragment). This only checks for one-way connection: the second
        fragment of self should overlap the first fragment of other."""
        both_frags = self.has_frags() and other.has_frags()
        both_signs = self.has_signs() and other.has_signs()
        compat_chroms = self.pos2.chrom == other.pos1.chrom
        if not (both_frags and both_signs):
            raise ValueError(
                "Cannot connect breakpoints without fragment or sign information."
            )
        compat_coords = False
        if self.pos2.sign != other.pos1.sign:
            if self.intersect(other) >= min_intersect:
                compat_coords = True

        return compat_chroms & compat_coords

    def intersect(self, other: BreakPoint) -> int:
        """Amount of overlap between two breakpoint's fragments"""
        if not (self.has_frags() and other.has_frags()):
            raise ValueError(
                "Cannot compute overlap if without fragment information."
            )
        # TODO: Could use signs to skip some comparisons
        intersect = 0
        for source in [self.frag1, self.frag2]:
            for target in [other.frag1, other.frag2]:
                intersect = max(intersect, source.intersect(target))
        return intersect
