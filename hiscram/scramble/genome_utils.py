"""
Utilities to generate work with genomes.
"""
import random
import numpy as np
import pandas as pd
import cooler
import pyfastx
from Bio import SeqIO, Seq
import json
from typing import Optional, Tuple
from dataclass import dataclass
from ..breakpoint import Fragment

@dataclass
class StruVar:
    """Represents a simple structural variant (DEL/INV/INS)."""
    chrom: str
    coord: int
    sv_type:  str


class GenomeMixer(object):
    """
    Handles genome edition through different types of structural variations.

    Examples
    --------
        mix = GenomeMixer("genome.fasta", "config.json", "profile="Dmel")
        mix.generate_sv()
        mix.edit_genome("new_genome.fasta")

    Attributes
    ----------
    fasta : pyfastx.Fasta
        Underlying genome.
    config : dict
        SV properties in a nested dictionary, loaded from a single profile in
        the config file.
    svs : List of StruVar
        Contains all structural variants, in the order they should be applied
        to restore the original genome.

    """

    def __init__(
        self,
        genome_path: str,
        config_path: str,
        config_profile: Optional[str] = None,
    ):
        self.config = self.load_profile(profile=config_profile)
        self._og_genome = pyfastx.Fasta(genome_path)
        self.svs = []

    @property
    def fasta(self):


    def add_sv(self, sv: StruVar):
        """
        Add a structural variant to the genome, and update coordinates of
        existing breakpoints accordingly
        """
        # Somehow implement complex SV by decomposing them into simple SVs:
        # TRA: INS [-> INV] -> DEL
        self.svs.append(sv)
        ...

    def load_profile(self, profile: Optional[str] = None):
        """
        Load SV profile from a JSON config file. The top level JSON object is a
        profile. A config file can have multiple profiles, but only one will be
        loaded. Each profile contains the 'SV_freq' property, which gives the
        number of SV per bp, and the `SV_types` object. The SV_types object
        contains one object per SV type. Each SV contains multiple properties.

        Parameters
        ----------
        profile : str
            Name of the profile to load in the config file.

        Returns
        -------
        dict :
            A dictionary of SV types structured like
            {"sv_type":{"property": value, ...}, ...}
        """
        config = json.load(open(self.config_path, "r"))
        if not profile:
            if len(config) > 1:
                print(
                    "You must specify a profile name if multiple profiles "
                    "appear in the JSON config file"
                )
                raise ValueError
            else:
                profile = config.keys()[0]

        return config[profile]


    def generate_sv(self) -> pd.DataFrame:
        """
        Generates random structural variations, based on the parameters loaded
        from the instance's config file.
        # NOTE: Currently only implemented for inversions and deletions.

        Returns
        -------
        pandas.DataFrame :
            A dataframe where each row is a SV. columns represent
            sv_type, chrom, start, end.
        """
        # Relative abundance of each event type (placeholder values)
        # rel_abun = {"INV": 8, "DEL": 400, "DUP": 60, "INS": 160, "CNV": 350}
        all_chroms_sv = []
        for chrom, size in self.chromsizes.items():
            n_sv = round(size * self.config["SV_freq"])
            chrom_sv = pd.DataFrame(np.empty((n_sv, 4)))
            chrom_sv.columns = ["sv_type", "chrom", "start", "end"]
            sv_count = 0
            for sv_name, sv_char in self.config["SV_types"].items():
                # multiply proportion of SV type by total SV freq desired to
                # get number of events of this type.
                n_event = round(n_sv * sv_char["prop"])

                # Safeguard rounding errors
                if sv_count + n_event > n_sv:
                    n_event -= (n_event + sv_count) - n_sv

                print("Generating {0} {1}".format(n_event, sv_name))
                for _ in range(n_event):
                    # Start position is random and length is picked from a normal
                    # distribution centered around mean length.
                    start = np.random.randint(size)
                    end = start + abs(
                        np.random.normal(
                            loc=sv_char["mean_size"], scale=sv_char["sd_size"]
                        )
                    )
                    # Make sure the inversion does not go beyond chromosome.
                    end = min(size, end)
                    chrom_sv.iloc[sv_count, :] = (sv_name, chrom, start, end)
                    sv_count += 1
            all_chroms_sv.append(chrom_sv)
        out_sv = pd.concat(all_chroms_sv, axis=0)
        out_sv.start, out_sv.end = (
            out_sv.start.astype(int),
            out_sv.end.astype(int),
        )
        # Randomize rows in SV table to mix the order of different SV types
        out_sv = out_sv.sample(frac=1).reset_index(drop=True)
        self.sv = out_sv

    def save_edited_genome(self, fasta_out: str):
        """
        Apply computed SVs to the sequence and store the edited sequence into
        the target file in fasta format.

        Coordinates in self.sv are updated as the genome is modified.

        Parameters
        ----------
        fasta_out : str
            Path where the edited genome will be written in fasta format.
        """
        with open(fasta_out, "w") as fa_out:
            for rec in SeqIO.parse(self.genome_path, format="fasta"):
                mutseq = Seq.MutableSeq(str(rec.seq))
                for row_num in range(self.sv.shape[0]):
                    row = self.sv.iloc[row_num, :]
                    sv_type = row.sv_type
                    chrom, start, end = row.chrom, int(row.start), int(row.end)
                    # NOTE: Only implemented for inversions for now.
                    if sv_type == "INV":
                        # Reverse complement to generate inversion
                        if chrom == rec.id:
                            mutseq[start:end] = Seq.reverse_complement(
                                mutseq[start:end]
                            )
                            # Update coordinates of other SVs in the INV region
                            mid = (end + start) // 2
                            starts = self.sv.eval(
                                "(chrom == @chrom) & (start >= @start) & (start <= @end)"
                            )
                            ends = self.sv.eval(
                                "(chrom == @chrom) & (end >= @start) & (end <= @end)"
                            )
                            self.sv.loc[starts, "start"] = (
                                mid + mid - self.sv.start[starts]
                            )
                            self.sv.loc[ends, "end"] = mid + mid - self.sv.end[ends]
                            # Make sure start is always lower than end
                            swap_mask = self.sv.start > self.sv.end
                            self.sv.loc[swap_mask, ["start", "end"]] = self.sv.loc[
                                swap_mask, ["end", "start"]
                            ].values
                    elif sv_type == "DEL":
                        if chrom == rec.id:
                            mutseq = mutseq[:start] + mutseq[end:]
                            # Shift coordinates on the right of DEL region
                            self.sv.loc[
                                (self.sv.chrom == chrom) & (self.sv.start >= start),
                                ["start", "end"],
                            ] -= (
                                end - start
                            )
                            self.sv.loc[self.sv.start < 0, "start"] = 0
                            self.sv.loc[self.sv.end < 0, "end"] = 0
                    else:
                        raise NotImplementedError("SV type not implemented yet.")
                self.sv.start = self.sv.start.astype(int)
                self.sv.end = self.sv.end.astype(int)
                # Discard SV that have been erased by others
                self.sv = self.sv.loc[(self.sv.end - self.sv.start) > 1, :]
                rec = SeqIO.SeqRecord(seq=mutseq, id=rec.id, description="")
                # Trim SV with coordinates > chrom size
                self.sv.loc[
                    (self.sv.chrom == chrom) & (self.sv.end >= len(mutseq)), "end"
                ] = (len(mutseq) - 1)
                SeqIO.write(rec, fa_out, format="fasta")


def pos_to_coord(
    clr: cooler.Cooler, sv_df: pd.DataFrame
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Converts start - end genomic positions from structural variations to breakpoints
    in matrix coordinates.

    Parameters
    ----------
    clr : cooler.Cooler
        The cooler object containing Hi-C data
    sv_df : pandas.DataFrame
        A dataframe containg the type and genomic start-end coordinates of
        strucural variations as given by generate_sv().

    Returns
    -------
    breakpoints : numpy.ndarray of int
        A N x 2 numpy array of numeric values representing X, Y coordinates of structural
        variations breakpoints in the matrix.
    labels : numpy.ndarray of str
        An N X 1 array of labels corresponding to SV type.
    """
    # Get coordinates to match binning
    res = clr.binsize
    sv_df.start = (sv_df.start // res) * res
    sv_df.end = (sv_df.end // res) * res
    # Put start and end in the same column, 1 row / breakpoint
    s_df = sv_df.loc[:, ["sv_type", "chrom", "start"]]
    s_df.rename(index=str, columns={"start": "pos"}, inplace=True)
    e_df = sv_df.loc[:, ["sv_type", "chrom", "end"]]
    e_df.rename(index=str, columns={"end": "pos"}, inplace=True)
    sv_df = pd.concat([s_df, e_df]).reset_index(drop=True)
    # Assign matrix coordinate (fragment index) to each breakpoint
    bins = clr.bins()[:]
    bins["coord"] = bins.index
    sv_frags = sv_df.merge(
        bins,
        left_on=["chrom", "pos"],
        right_on=["chrom", "start"],
        how="left",
    )
    breakpoints = np.vstack([sv_frags.coord, sv_frags.coord]).T
    breakpoints.astype(int)
    labels = np.array(sv_frags.sv_type.tolist())
    return breakpoints, labels


def get_random_region(
        fasta: pyfastx.Fasta,
        region_size: int = 1000
        ) -> Fragment:
    """Return a random region from the genome (chrom, start, end)."""
    # Exclude chromosomes smaller than slice_size
    chrom_sizes = {
        seq.name: len(seq) for seq in fasta if len(seq) < slice_size
    }

    # Get list of valid chromosomes
    chrom_names = list(chrom_sizes.keys())

    # Pick a random region of slice_size bp in a random chromosome and write it
    picked_chrom = random.choices(
        chrom_sizes.keys(),
        weights=chrom_sizes.values(),
        k=1,
    )[0]
    start_slice = int(
        np.random.randint(low=0, high=chrom_sizes[picked_chrom] - slice_size, size=1)
    )
    end_slice = int(start_slice + slice_size)

    return Fragment(picked_chrom, start_slice, end_slice)


def save_genome_slice(
        fasta: pyfastx.Fasta,
        region: Fragment,
        out_path: str,
    ):
    """
    Given an input fasta file, slice a random region of a random chromosome and
    save it into a new fasta file.

    Parameters
    ----------
    fasta:
        Genome from which to take the slice.
    region:
        Coordinates of the region to extract, in basepairs.
    out_path:
        Coordinates of the region to extract.

    """

    with open(out_path, "w") as sub_handle:
        seq = fasta[region.chrom]
        sliced_seq = seq[region.start:region.end]
        sub_handle.write(f'>{region[0]}\n{sliced_seq}\n')


