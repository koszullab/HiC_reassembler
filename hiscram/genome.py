"""
Datastructures representing a genome as a collection of ragments linked by breakpoints.
A number of methods allow to introduce structural variation in the genome, and retrieve
the altered sequence.
"""

from typing import List, Dict, Iterator
import pyfastx
import numpy as np
from hiscram.regions import BreakPoint, Fragment, Position


class Chromosome:
    """Representation of a chromosome as a collection of fragments.
    Each fragment represents a (0-based, right-open) region of the original genome."""

    def __init__(self, name: str, length: int):
        self.name = name
        self.frags = [Fragment(self.name, 0, length)]
        self.breakpoints = []

    def __len__(self):
        """Returns the total chromosome length."""
        return sum(len(frag) for frag in self.frags)

    @property
    def boundaries(self):
        """Get array of fragment boundaries, from the start to the end
        of the chromosome."""
        positions = np.cumsum([0] + [len(frag) for frag in self.frags])
        return positions

    def clean_frags(self):
        """Purge 0-length fragments."""
        self.frags = [frag for frag in self.frags if len(frag)]

    def insert(self, position: int, frag_ins: Fragment):
        """Updates fragments by inserting a sequence in the chromosome."""
        bounds = self.boundaries
        if position in bounds:
            # Insertion right between two fragments, add a fragment.
            frag_id = np.where(bounds == position)[0][0]
            self.frags.insert(frag_id, frag_ins)
        else:
            # Insertion inside a fragment, split it and add fragment in between.
            frag_id = max(np.searchsorted(bounds, position) - 1, 0)
            frag_l, frag_r = self.frags.pop(frag_id).split(position - bounds[frag_id])
            for frag in [frag_r, frag_ins, frag_l]:
                self.frags.insert(frag_id, frag)

    def invert(self, start: int, end: int):
        """Updates fragments by inverting a portion of the chromosome."""
        bounds = self.boundaries
        frag_start = max(np.searchsorted(bounds, start) - 1, 0)
        frag_end = max(np.searchsorted(bounds, end) - 1, 0)
        start_dist = start - bounds[frag_start]
        end_dist = end - bounds[frag_end]

        # Inversion inside a single frag.: Split it in 3 and invert middle.
        if frag_end == frag_start:
            inv_size = end - start
            frag_l, frag_mr = self.frags.pop(frag_start).split(start_dist)
            frag_m, frag_r = frag_mr.split(inv_size)
            frag_m.flip()
            for frag in [frag_r, frag_m, frag_l]:
                self.frags.insert(frag_start, frag)
        else:
            # Split fragment where inversion starts and flip right part.
            start_l, start_r = self.frags.pop(frag_start).split(start_dist)
            start_r.flip()
            for frag in [start_r, start_l]:
                self.frags.insert(frag_end, frag)
            # If fragments are entirely in the inversion, invert and flip them.
            for frag_id in range(frag_start + 1, frag_end):
                inv_frag = self.frags.pop(frag_id)
                inv_frag.flip()
                self.frags.insert(frag_id)

            # Split fragment where inversion ends and flip left part.
            end_l, end_r = self.frags.pop(frag_end).split(end_dist)
            end_l.flip()
            for frag in [end_r, end_l]:
                self.frags.insert(frag_end, frag)
        self.clean_frags()

    def delete(self, start: int, end: int):
        """Updates fragments by deleting a portion of the chromosome."""
        bounds = self.boundaries
        frag_start = max(np.searchsorted(bounds, start) - 1, 0)
        frag_end = max(np.searchsorted(bounds, end) - 1, 0)
        del_size = end - start
        start_dist = start - self.frags[frag_start].start
        # Deletion contained in a single fragment: split it and trim right part
        if frag_end == frag_start:
            start_l, start_r = self.frags.pop(frag_start).split(start_dist)
            start_r.start += del_size
            for frag in [start_r, start_l]:
                self.frags.insert(frag_start, frag)
        # Deletion spans multiple fragments
        else:
            # Deletion starts in frag, end gets trimmed
            self.frags[frag_start].end = start

            # Fragments contained in deletion disappear
            for frag_id in range(frag_start + 1, frag_end):
                if self.frags[frag_id].end < end:
                    self.frags[frag_id].end = start
                if self.frags[frag_id].start < end:
                    self.frags[frag_id].start = start

            # Deletion ends in frag, trim left side
            self.frags[frag_end].start = end
        self.clean_frags()

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
            p1 = Position(self.name, frag1.end, not frag1.is_reverse)
            p2 = Position(self.name, frag2.start, frag2.is_reverse)
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
            self.chroms[target_chrom].invert(target_pos, target_pos + frag_size)
        self.chroms[source_region.chrom].delete(source_region.start, source_region.end)

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
