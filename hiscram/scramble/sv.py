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
from hiscram.genome import Genome
from hiscram.breakpoint import Fragment


class Mixer(object):
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
        self.genome = Genome(pyfastx.Fasta(genome_path))
        self.svs = []

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

        Returns
        -------
        pandas.DataFrame :
            A dataframe where each row is a SV. columns represent
            sv_type, chrom, start, end.
        """
        # Relative abundance of each event type (placeholder values)
        # rel_abun = {"INV": 8, "DEL": 400, "DUP": 60, "INS": 160, "CNV": 350}
        all_chroms_sv = []
        for chrom, size in self.genome.chromsizes.items():
            n_sv = round(size * self.config["SV_freq"])
            sv_count = 0
            # TODO: Randomize order of SVs (currently all INV, then all DEL...)
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
                    size = abs(
                        np.random.normal(
                            loc=sv_char["mean_size"], scale=sv_char["sd_size"]
                        )
                    )
                    region = get_random_region(self.genome, size)
                    target_pos = get_random_region(self.genome, 1)
                    # Make sure the inversion does not go beyond chromosome.
                    end = min(size, end)
                    if sv_name == "DEL":
                        self.genome.delete(
                            region.chrom, region.start, region.end
                        )
                    elif sv_name == "INS":
                        self.genome.insert(
                            target_pos.chrom, target_pos.start, region
                        )
                    elif sv_name == "INV":
                        self.genome.invert(
                            region.chrom, region.start, region.end
                        )
                    elif sv_name == "TRA":
                        inv = np.random.choice([False, True])
                        self.genome.translocate(
                            target_pos.chrom, target_pos.start, region, inv
                        )

                    sv_count += 1

    def get_breakpoints(self):
        ...


def get_random_region(genome: Genome, region_size: int = 1000) -> Fragment:
    """Return a random region from the genome (chrom, start, end)."""
    # Exclude chromosomes smaller than slice_size
    chrom_sizes = {
        name: size
        for name, size in genome.chromsizes.items()
        if size > region_size
    }

    # Pick a random region of slice_size bp in a random chromosome and write it
    picked_chrom = random.choices(
        chrom_sizes.keys(),
        weights=chrom_sizes.values(),
        k=1,
    )[0]
    start_slice = int(
        np.random.randint(
            low=0, high=chrom_sizes[picked_chrom] - region_size, size=1
        )
    )
    end_slice = int(start_slice + region_size)

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
        sliced_seq = seq[region.start : region.end]
        sub_handle.write(f">{region[0]}\n{sliced_seq}\n")
