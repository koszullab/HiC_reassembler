"""
Datastructures representing a genome as a collection of ragments linked by breakpoints.
A number of methods allow to introduce structural variation in the genome, and retrieve
the altered sequence.
"""
from ast import Break
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
        coordinate falls. Return format is (id, (start, end))."""
        bounds = self.boundaries
        if coord >= bounds[-1]:
            raise ValueError(f"Coordinate out of bounds: {self.name}:{coord}")
        frag_id = max(0, np.searchsorted(bounds, coord, side="right") - 1)
        return (frag_id, (bounds[frag_id], bounds[frag_id + 1]))
    
    def get_frags_between(self, start: int, end:int) -> List[Fragment]:
        """Returns a list of the fragments between 2 positions in the chromosome"""
        result = []
        bounds = self.boundaries
        if start < 0 or end > bounds[-1]:
            raise ValueError(f"Coordinate out of bounds: {self.name}:{start}-{end}")

        frag_id_left, (start_bound_left, end_bound_left) = self.get_frag_bounds(start)
        frag_id_right, (start_bound_right, end_bound_right) = self.get_frag_bounds(end)

        # If the start and end falls in the same fragment
        if(frag_id_left == frag_id_right):
            return [Fragment(
                self.name,
                self.frags[frag_id_left].start + (start - start_bound_left),
                self.frags[frag_id_left].start + (end - start_bound_left),
                self.frags[frag_id_left].is_reverse)]
        # Else
        # Start falls in a fragment, trim it
        print("start index: ", start)
        print("start bound: ", self.get_frag_bounds(start))
        left_frag = self.frags[frag_id_left]
        print("left_frag = ", left_frag)

        left_frag_trimmed = Fragment(
            self.name,
            left_frag.start + (start - start_bound_left),
            left_frag.end,
            left_frag.is_reverse
        )
        print("left frag trimmed = ", left_frag_trimmed)
        result.append(left_frag_trimmed)

        # Ends falls in a fragment, trim it
        right_frag = self.frags[frag_id_right]
        print("right_frag = ", right_frag)

        right_frag_trimmed = Fragment(
            self.name,
            right_frag.start,
            right_frag.start + end - start_bound_right,
            right_frag.is_reverse
        )
        print("right frag trimmed = ", right_frag_trimmed)
        # Add all fragments between the first and the last
        for frag_i in range(frag_id_left + 1, frag_id_right):
            frag = Fragment(
                self.name,
                self.frags[frag_i].start,
                self.frags[frag_i].end,
                self.frags[frag_i].is_reverse
            )
            result.append(frag)
        
        # Finally add last fragment
        result.append(right_frag_trimmed)

        return result
    
    def insert_frags(self, position: int, frags: List[Fragment]):
        insert_pos = position
        for i in range(len(frags)):
            self.insert(insert_pos, frags[i])
            insert_pos += len(frags[i])

    def clean_frags(self):
        """Purge 0-length fragments."""
        self.frags = [frag for frag in self.frags if len(frag)]

    def insert(self, position: int, frag_ins: Fragment):
        """Updates fragments by inserting a sequence in the chromosome."""
        bounds = self.boundaries
        # Append after the end of chromosome
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
        bp = BreakPoint(Position(self.name, frag_ins.start,),Position(self.name, frag_ins.end), "INS")
        self.breakpoints.append(bp)
        

    def invert(self, start: int, end: int):
        """Updates fragments by inverting a portion of the chromosome.
        The interval is 0-based and right open [start;end[."""
        s_frag_id, (s_frag_start, _) = self.get_frag_bounds(start)
        e_frag_id, (e_frag_start, _) = self.get_frag_bounds(end - 1)
        s_start_dist = start - s_frag_start
        e_start_dist = end - e_frag_start

        # Inversion inside a single frag.: Split it in 3 and invert middle.
        if s_frag_id == e_frag_id:
            inv_size = end - start
            frag_l, frag_mr = self.frags.pop(s_frag_id).split(s_start_dist)
            frag_m, frag_r = frag_mr.split(inv_size)
            frag_m.flip()
            for frag in [frag_r, frag_m, frag_l]:
                self.frags.insert(s_frag_id, frag)
        else:
            # Split fragment where inversion starts, we'll flip the right part.
            start_l, start_r = self.frags.pop(s_frag_id).split(s_start_dist)
            for frag in [start_r, start_l]:
                self.frags.insert(s_frag_id, frag)
            s_frag_id += 1
            e_frag_id += 1
            # Split fragment where inversion ends we'll flip the left part.
            end_l, end_r = self.frags.pop(e_frag_id).split(e_start_dist)
            for frag in [end_r, end_l]:
                self.frags.insert(e_frag_id, frag)
            e_frag_id += 1
            # If fragments are contained in the inversion, invert and flip them.
            for frag_id in range(s_frag_id, e_frag_id):
                self.frags[frag_id].flip()
            self.frags[s_frag_id:e_frag_id] = self.frags[
                e_frag_id - 1 : s_frag_id - 1 : -1
            ]
        self.clean_frags()
        bp = BreakPoint(Position(self.name,start),Position(self.name,end), "INV")
        self.breakpoints.append(bp) 

    def delete(self, start: int, end: int):
        """Updates fragments by deleting a portion of the chromosome.
        The interval is 0-based and right open [start;end[."""
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
                self.frags[frag_id].end = curr_start

            from copy import copy

            ori_end = copy(self.frags[e_frag_id])
            # Deletion ends in frag, trim left side
            self.frags[e_frag_id].start = self.frags[e_frag_id].end - end_dist

        bp = BreakPoint(Position(self.name,start),Position(self.name,end), "DEL")
        self.breakpoints.append(bp)

    def get_seq(self, fasta: pyfastx.Fasta) -> Iterator[str]:
        """Retrieve the chromosome sequence, as a generator yielding
        the sequence by fragment."""
        self.clean_frags()
        for frag in self.frags:
            strand = "-" if frag.is_reverse else "+"
            # Note: fasta.fetch is 1-based...
            yield fasta.fetch(
                frag.chrom,
                (int(frag.start + 1), (frag.end)),
                strand=strand,
            )

    def get_breakpoints(self) -> Iterator[BreakPoint]:
        """Retrieve a generator yielding breakpoints in the chromosome."""
        self.clean_frags()
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
    
    def translocate(
        self,
        target_chrom: str,
        target_pos: int,
        source_chrom: str,
        source_start: int,
        source_end:int,
        invert: bool = False,
    ):
        """Move a portion of the chromosome to another genomic position.
        Is not limited at 1 fragment of original genome, like the previous function."""
        frags_to_insert = self.chroms[source_chrom].get_frags_between(source_start, source_end)
        insert_size = source_end - source_start
        self.chroms[target_chrom].insert_frags(target_pos, frags_to_insert)
        if(invert):
            self.chroms[target_chrom].invert(target_pos, target_pos + insert_size)
        if(target_chrom == source_chrom):
            if(target_pos > source_end):
                self.chroms[source_chrom].delete(source_start, source_end)
            else:
                self.chroms[source_chrom].delete(source_start + insert_size, source_end + insert_size)
        else:
            self.chroms[source_chrom].delete(source_start, source_end)

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