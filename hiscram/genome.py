"""
Datastructures representing a genome as a collection of ragments linked by breakpoints.
A number of methods allow to introduce structural variation in the genome, and retrieve
the altered sequence.
"""

from typing import List, Dict, Iterator, Tuple
import pyfastx
import numpy as np
from hiscram.regions import BreakPoint, Fragment, Position


class Chromosome:
    """Representation of a chromosome as a collection of fragments.
    Each fragment represents a (0-based, right-open) region of the original genome."""

    def __init__(self, name: str, length: int):
        self.name = name
        # list rather than np.array, due to better insert/append performance
        self.frags = [Fragment(self.name, 0, length)]
        self.breakpoints = []

    def __len__(self):
        """Returns the total chromosome length."""
        return self.boundaries[-1]

    @property
    def boundaries(self):
        """Get array of fragment boundaries, from the start to the end
        of the chromosome."""
        # Memorize whether fragments have changed to avoid recomputing the same
        # values.
        frags_hash = hash(tuple(self.frags))
        try:
            if self._frags_hash == frags_hash:
                changed = False
            else:
                changed = True
        # On first access, required attrs are generated
        except AttributeError:
            self._frags_hash = frags_hash
            changed = True
        if changed:
            self._bds = np.cumsum([0] + [len(frag) for frag in self.frags])
        return self._bds

    def get_frag_bounds(self, coord: int) -> Tuple[int, Tuple[int, int]]:
        """Returns the index and boundaries of the fragment in which input
        coordinate falls. Return format is (id, start, end)"""
        bounds = self.boundaries
        if coord >= bounds[-1]:
            raise ValueError(f"Coordinate out of bounds: {self.name}:{coord}")
        frag_id = max(0, np.searchsorted(bounds, coord, side="right") - 1)
        return (frag_id, (bounds[frag_id], bounds[frag_id + 1]))

    def clean_frags(self):
        """Purge 0-length fragments."""
        self.frags = [frag for frag in self.frags if len(frag)]

    def insert(self, position: int, frag_ins: Fragment):
        """Updates fragments by inserting a sequence in the chromosome."""
        bounds = self.boundaries
        if position == len(self):
            frag_id = len(bounds)
        else:
            frag_id, (frag_start, _) = self.get_frag_bounds(position)
        if position in bounds:
            # Insertion right between two fragments, add a fragment.
            self.frags.insert(frag_id, frag_ins)
        else:
            # Insertion inside a fragment, split it and add fragment in between.
            frag_l, frag_r = self.frags.pop(frag_id).split(
                position - frag_start
            )
            for frag in [frag_r, frag_ins, frag_l]:
                self.frags.insert(frag_id, frag)

    def invert(self, start: int, end: int):
        """Updates fragments by inverting a portion of the chromosome."""
        s_frag_id, (s_frag_start, _) = self.get_frag_bounds(start)
        e_frag_id, (_, e_frag_end) = self.get_frag_bounds(end - 1)
        start_dist = start - s_frag_start
        end_dist = e_frag_end - end

        # Inversion inside a single frag.: Split it in 3 and invert middle.
        if s_frag_id == e_frag_id:
            inv_size = end - start
            frag_l, frag_mr = self.frags.pop(s_frag_id).split(start_dist)
            frag_m, frag_r = frag_mr.split(inv_size)
            frag_m.flip()
            for frag in [frag_r, frag_m, frag_l]:
                self.frags.insert(s_frag_id, frag)
        else:
            # Split fragment where inversion starts and flip right part.
            start_l, start_r = self.frags.pop(s_frag_id).split(start_dist)
            start_r.flip()
            for frag in [start_r, start_l]:
                self.frags.insert(e_frag_id, frag)
            # If fragments are entirely in the inversion, invert and flip them.
            for frag_id in range(s_frag_id + 1, e_frag_id):
                inv_frag = self.frags.pop(frag_id)
                inv_frag.flip()
                self.frags.insert(frag_id)

            # Split fragment where inversion ends and flip left part.
            end_l, end_r = self.frags.pop(e_frag_id).split(end_dist)
            end_l.flip()
            for frag in [end_r, end_l]:
                self.frags.insert(e_frag_id, frag)
        self.clean_frags()

    def delete(self, start: int, end: int):
        """Updates fragments by deleting a portion of the chromosome."""
        s_frag_id, (s_frag_start, _) = self.get_frag_bounds(start)
        e_frag_id, (_, e_frag_end) = self.get_frag_bounds(end - 1)
        del_size = end - start
        start_dist = start - s_frag_start
        end_dist = e_frag_end - end
        # Deletion contained in a single fragment: split it and trim right part
        if e_frag_id == s_frag_id:
            start_l, start_r = self.frags.pop(s_frag_id).split(start_dist)
            start_r.start += del_size
            for frag in [start_r, start_l]:
                self.frags.insert(s_frag_id, frag)
        # Deletion spans multiple fragments
        else:
            # Deletion starts in frag, end gets trimmed
            self.frags[s_frag_id].end = (
                self.frags[s_frag_id].start + start_dist
            )

            # Fragments contained in deletion disappear
            for frag_id in range(s_frag_id + 1, e_frag_id):
                curr_start = self.frags[frag_id].start
                if self.frags[frag_id].end < end:
                    self.frags[frag_id].end = curr_start
                if self.frags[frag_id].start < end:
                    self.frags[frag_id].start = curr_start

            from copy import copy

            ori_end = copy(self.frags[e_frag_id])
            # Deletion ends in frag, trim left side
            self.frags[e_frag_id].start = self.frags[e_frag_id].end - end_dist

    def get_seq(self, fasta: pyfastx.Fasta) -> Iterator[str]:
        """Retrieve the chromosome sequence, as a generator yielding
        the sequence by fragment."""
        for frag in self.frags:
            strand = "-" if frag.is_reverse else "+"
            # Note: fasta.fetch is 1-based...
            yield fasta.fetch(
                frag.chrom, (int(frag.start + 1), (frag.end)), strand=strand
            )

    def get_breakpoints(self) -> Iterator[BreakPoint]:
        """Retrieve a generator yielding breakpoints in the chromosome."""
        for frag_id in range(0, len(self.frags) - 1):
            frag1, frag2 = self.frags[frag_id : frag_id + 2]
            p1 = Position(frag1.chrom, frag1.end, not frag1.is_reverse)
            p2 = Position(frag2.chrom, frag2.start, frag2.is_reverse)
            breakpoint = BreakPoint(p1, p2)
            breakpoint.frag1 = frag1
            breakpoint.frag2 = frag2
            yield breakpoint


class Genome:
    """Collection of chromosomes allowing complex SVs such as translocations."""

    def __init__(self, fasta: pyfastx.Fasta):
        self.fasta = fasta
        self.chroms = {}
        for seq in fasta:
            self.chroms[seq.name] = Chromosome(seq.name, len(seq))

    def __len__(self):
        return sum([len(chrom) for chrom in self.chroms.values()])

    @property
    def chromsizes(self) -> Dict[str, int]:
        chromsizes = {}
        for chrom in self.chroms.values():
            chromsizes[chrom.name] = len(chrom)
        return chromsizes

    def delete(self, chrom: str, start: int, end: int):
        """Delete a genomic segment."""
        self.chroms[chrom].delete(start, end)

    def insert(self, chrom: str, position: int, frag: Fragment):
        """Insert a new genomic segment."""
        self.chroms[chrom].insert(position, frag)

    def invert(self, chrom: str, start: int, end: int):
        """Invert (i.e. flip) a genomic segment."""
        self.chroms[chrom].invert(start, end)

    def translocate(
        self,
        target_chrom: str,
        target_pos: int,
        source_region: Fragment,
        invert: bool = False,
    ):
        """Move a genomic segment to another genomic position."""
        frag_size = source_region.end - source_region.start
        self.chroms[target_chrom].insert(target_pos, source_region)
        if invert:
            self.chroms[target_chrom].invert(
                target_pos, target_pos + frag_size
            )
        self.chroms[source_region.chrom].delete(
            source_region.start, source_region.end
        )

    def duplicate(
        self,
        target_chrom: str,
        target_pos: int,
        source_region: Fragment,
        invert: bool = False,
    ):
        """Copy a genomic segment to another genomic position."""
        frag_size = source_region.end - source_region.start
        self.chroms[target_chrom].insert(target_pos, source_region)
        if invert:
            self.chroms[target_chrom].invert(
                target_pos, target_pos + frag_size
            )

    def get_seq(self) -> Dict[str, Iterator[str]]:
        """Retrieve the genomic sequence of each chromosome. Each chromosome's
        sequence is returned as a generator of (lazily-retrieved) fragment sequences."""
        seqs = {}
        for chrom in self.chroms.values():
            seqs[chrom.name] = chrom.get_seq(self.fasta)
        return seqs

    def get_breakpoints(self) -> List[BreakPoint]:
        """Retrieve the list of all breakpoints between genomic fragments"""
        bps = []
        for chrom in self.chroms.values():
            bps += [br for br in chrom.get_breakpoints()]
        return bps
