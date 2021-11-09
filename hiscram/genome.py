from typing import Optional
from dataclasses import dataclass
import numpy as np
import hiscram.breakpoint as bp


class Chromosome:
    """Representation of a chromosome as a collection of breakpoints and fragments"""

    def __init__(self, name: str, length: int):
        self.name = name
        self.frags = [bp.Fragment(self.name, 0, length)]
        self.breakpoints = []

    def __len__(self):
        """Returns the total chromosome length."""
        return sum(len(frag) for frag in self.frags)

    @property
    def starts(self):
        """Get array of start positions from all fragments."""
        starts = [frag.start for frag in self.frags]
        return starts

    @property
    def ends(self):
        """Get array of end positions from all fragments."""
        ends = [frag.end for frag in self.frags]
        return ends

    def insert(self, position: int, size: int, seq: Optional[str]=None):
        frag_ins = np.searchsorted(self.starts, position)
        for frag_id in range(frag_ins+1)

    def invert(self, start: int, end: int):
        frag_start = np.searchsorted(self.starts, start)
        frag_end = np.searchsorted(self.ends, end)

        # Fragment where deletion starts is split and right part is flipped
        start_l, start_r = self.frags.pop(frag_start).split(start - frag_start.start)
        self.frags.insert(start_l, frag_start)
        self.frags.insert(start_r.flip(), frag_start + 1)

        # Fragments inside deletions are flipped
        for frag_id in range(frag_start + 1, frag_end):
            self.frags[frag_id].flip()
        
        # Split last fragment where inversion starts and flip left part
        end_l, end_r = self.frags.pop(frag_end).split(frag_end - end)
        self.frags.insert(end_l.flip(), frag_end)
        self.frags.insert(end_r, frag_end + 1)

    def delete(self, start: int, end: int):
        frag_start = np.searchsorted(self.starts, start)
        frag_end = np.searchsorted(self.ends, end)
        del_size = end - start
        # Deletion starts in frag, end gets trimmed
        self.frags[frag_start].end = start
        # Fragments contained in deletion disappear
        for frag_id in range(frag_start + 1, frag_end):
            self.frags[frag_id].start = start
            self.frags[frag_id].end = start
        # Deletion ends in frag, different shifts
        self.frags[frag_end].start -= start
        self.frags[frag_end].end -= del_size
        # Fragments downstream of deletion -> shifted
        for frag in self.frags[frag_end+1:]:
            frag.start -= del_size
            frag.end -= del_size
        # Update fragments to delete those that were removed
        self.frags = [frag for frag in self.frags if len(frag)]


class Genome:
    """Collection of chromosomes allowing complex SVs like translocations."""
    def __init__(self, fasta: pyfastx.Fasta):
        self.fasta = fasta
        self.chroms = {}
        for seq in fasta:
            self.chroms[seq.name] = Chromosome(seq.name, len(seq))
    
    def delete(self, chrom: str, start: int, end: int):
        self.chroms[chrom].delete(start, end)
    
    def insert(self, chrom: str, position: int, size: int, seq: Optional[str]=None):
        self.chroms[chrom].insert(position, size, seq)
    
    def invert(self, chrom: str, start: int, end: int):
        self.chroms[chrom].invert(start, end)

    def translocate(self, region: bp.Fragment, target: bp.Position, invert: bool = False):
        frag_size = region.end - region.start
        self.chroms[target.chrom].insert(target.coord, frag_size)
        if invert:
            self.chroms[target.chrom].invert(target.coord, target.coord + frag_size)
        self.chroms[region.chrom].delete(region.start, region.end)
    