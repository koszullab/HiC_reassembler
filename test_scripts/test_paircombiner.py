from email.policy import default
from tempfile import tempdir
from hiscram.detector.pair_combiner import PairCombiner
import click
import pathlib
import numpy as np
from os import path, mkdir
import sys
from matplotlib import pyplot as plt
from hiscram.genome import Genome
import pyfastx

@click.command()
@click.option(
    "--number",
    "-n",
    type=int,
    default=1,
)
@click.option(
    "--chrom",
    "-c",
    type=str,
    required=True,
)
@click.option(
    "--binsize",
    "-b",
    type=int,
    default=1000,
    help="Size of the bins used for the matrices"
)
@click.option(
    "--readlength",
    "-rl",
    type=int,
    default=36,
)
@click.option(
    "--tmpdir",
    "-T",
    type=click.Path(path_type=pathlib.Path),
    default="./tmp"
)
@click.argument("datadir", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("name", type=str, default="scrambled")
@click.argument("index", type=int, default=None)
def test_pair_combiner(number,chrom, binsize, datadir, name,readlength,tmpdir,index):
    """
    Tests the association of known breakpoints based on the pairfile.
    """
    scrambledir = datadir / f"scrambling"
    try:
        mkdir(tmpdir)
    except:
        pass

    # Retrieve the number of scrambling matrices available
    n_dir = 0
    while path.isdir(scrambledir / f"{name}_{n_dir}"):
        n_dir += 1
    
    for i in range(number):
        # Select a random matrix
        print(i)
        breakpoints = []
        while True:
            if index is None:
                r = int(np.random.uniform(low=0, high=n_dir))
            else:
                r = index
            cool_mat = str(scrambledir/ f"{name}_{r}" / f"{name}_{r}.cool")
            # Retreive the breakpoints
            with open(scrambledir / f"{name}_{r}" / f"breakpoints_{r}.txt") as bp:
                for line in bp:
                    # L = line.strip().split()
                    sv_type = ""
                    # if sv_type not in tested_SVs : raise ValueError()
                    # line_breakpoints = L[1].split(':')[1].split('-')
                    breakpoints.append(int(line))
                    # breakpoints += [int(int(b)/binsize) for b in line_breakpoints]
            break
        breakpoints = np.unique(breakpoints)

        paircombiner = PairCombiner(
            chrom, 
            breakpoints,
            scrambledir / f"{name}_{r}" / "tmp" / f"{name}_{r}.valid_idx.pairs",
            binsize,
            readlength,
            str(tmpdir),
            )
        paircombiner.create_score_matrix()
        print(f"breakpoints_{r} -> {paircombiner.bp_to_associate}")
        paircombiner.solve_matrix(debug=True)
        paircombiner.print_sv()

        
        genome = Genome(pyfastx.Fasta(f"{scrambledir}/{name}_{r}/mod_genome_{r}.fa"))
        for frag in genome.chroms[chrom].frags:
            print(frag)


if __name__ == "__main__":
    test_pair_combiner()