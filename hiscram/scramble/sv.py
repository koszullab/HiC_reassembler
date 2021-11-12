"""
Utilities to generate work with genomes.
"""
import sys
import random
import json
from typing import Dict, Union, Any
import numpy as np
import pyfastx
from typing import Optional, Tuple
from hiscram import SCRAMBLE_CONFIG_PATH
from hiscram.genome import Genome
from hiscram.regions import Fragment


class Mixer(object):
    """
    Handles genome edition through different types of structural variations.

    Examples
    --------
        mix = Mixer("genome.fasta", "config.json", "profile="Dmel")
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
        config_path: str = SCRAMBLE_CONFIG_PATH,
        config_profile: Optional[str] = None,
    ):
        self.config = self.load_profile(config_path, profile=config_profile)
        self.genome = Genome(pyfastx.Fasta(genome_path))
        self.svs = []

    def load_profile(self, path: str, profile: Optional[str] = None) -> Dict[str, Any]:
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
        config = json.load(open(path, "r"))
        if not profile:
            if len(config) > 1:
                print(
                    "You must specify a profile name if multiple profiles "
                    "appear in the JSON config file"
                )
                raise ValueError(
                    f"No profile specified. Valid profiles are: {', '.join(config.keys())}"
                )
            else:
                profile = config.keys()[0]

        return config[profile]

    def generate_sv(self):
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
        n_sv = round(len(self.genome) * self.config["SV_freq"])
        sv_count = 0
        # multiply proportion of each SV type by total desired SV freq to
        # get total number of events of each type.
        n_events = {
            name: np.floor(n_sv * sv["prop"])
            for name, sv in self.config["SV_types"].items()
        }

        print(
            f"Generating {', '.join([str(num) + ' ' + sv for sv, num in n_events.items()])}"
        )
        for event_type in np.random.choice(n_events.keys(), n_events.values()):
            event_cfg = self.config["SV_types"][event_type]
            # Start position is random and length is picked from a normal
            # distribution centered around mean length.
            size = abs(
                np.random.normal(loc=event_cfg["mean_size"], scale=event_cfg["sd_size"])
            )
            region = get_random_region(self.genome, size)
            target_pos = get_random_region(self.genome, 1)
            # Make sure the inversion does not go beyond chromosome.
            end = min(size, end)
            self.svs.append((event_type, region, target_pos))
            if event_type == "DEL":
                self.genome.delete(region.chrom, region.start, region.end)
            elif event_type == "INS":
                self.genome.insert(target_pos.chrom, target_pos.start, region)
            elif event_type == "INV":
                self.genome.invert(region.chrom, region.start, region.end)
            elif event_type == "TRA":
                inv = np.random.choice([False, True])
                self.genome.translocate(target_pos.chrom, target_pos.start, region, inv)
            else:
                raise NotImplementedError("Unknown SV type.")

            sv_count += 1


def get_random_region(genome: Genome, region_size: int = 1000) -> Fragment:
    """Return a random region from the genome (chrom, start, end)."""
    # Exclude chromosomes smaller than slice_size
    chrom_sizes = {
        name: size for name, size in genome.chromsizes.items() if size > region_size
    }

    # Pick a random region of slice_size bp in a random chromosome and write it
    picked_chrom = random.choices(
        chrom_sizes.keys(), weights=chrom_sizes.values(), k=1,
    )[0]
    start_slice = int(
        np.random.randint(low=0, high=chrom_sizes[picked_chrom] - region_size, size=1)
    )
    end_slice = int(start_slice + region_size)

    return Fragment(picked_chrom, start_slice, end_slice)


def save_genome_slice(
    fasta: pyfastx.Fasta, region: Fragment, out_path: str,
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
