from distutils import core
from multiprocessing.sharedctypes import Value
from tkinter.tix import DirSelectDialog
from hiscram.detector.bamdetectormax import Bamdetectormax
from hiscram.detector.bam_functions import (
    bam_region_coverage, bam_region_read_ends, check_gen_sort_index, parse_ucsc_region)

import click
from matplotlib import pyplot as plt
import numpy as np
import pathlib
import pysam as ps
import typing
from os import path
import pyfastx

from scripts.detect import detect

@click.command()
@click.option(
    "--chrom",
    "-c",
    type=str,
    required=True,
)
@click.option(
    "--verbatim",
    "-v",
    is_flag=True,
    help="Prints details on each test",
)
@click.option(
    "--name",
    "-n",
    type=str,
    default="scrambled",
)
@click.argument("datadir", type=click.Path(exists=True,path_type=pathlib.Path))
def test_bam_policy(datadir,name, chrom,verbatim):
    """
    Tests the detection of breakpoints at base pair precision. Plots the percentage of True Positives and the median miss.
    """

    scrambledir = datadir / f"scrambling"
    winsize = 4000

    policy_detected_bps = 0
    total_bps = 0
    deviations = []

    # Retrieve the number of scrambling matrices available
    n_dir = 0
    while path.isdir(scrambledir / f"{name}_{n_dir}"):
        n_dir += 1
    
    for i in range(n_dir):
        breakpoints = []
        # Retreive the breakpoints
        with open(scrambledir / f"{name}_{i}" / f"breakpoints_{i}.txt") as bp:
            for line in bp:
                L = line.strip()
                breakpoints += [int(L)]
        breakpoints = np.unique(breakpoints)
        total_bps += len(breakpoints)

        # Load the alignment files
        try:
            bamfile_for = ps.AlignmentFile(scrambledir / f"{name}_{i}" / "tmp" / f"{name}_{i}.for.bam")
            bamfile_rev = ps.AlignmentFile(scrambledir / f"{name}_{i}" / "tmp" / f"{name}_{i}.rev.bam")
            bam_sorted_for = check_gen_sort_index(bamfile_for)
            bam_sorted_rev = check_gen_sort_index(bamfile_rev)
        except:
            continue
        
        # Load genome
        genome = str(scrambledir / f"{name}_{i}" / f"mod_genome_{i}.fa")


        # For each bp, apply policy on a window around this bp
        maxindexes = []
        detector = Bamdetectormax(chrom,1000,bam_sorted_for, bam_sorted_rev,genome)

        for bp in breakpoints:
            start = bp - winsize//2
            end = bp + winsize//2
            position = detector.old_predict(start, end, divided_by_coverage=True)
            maxindexes.append(position)
            

        # Record results
        positive_trigger_length = 5
        file_deviation = []

        for bp in breakpoints:
            deviances = [abs(bp - prediction) for prediction in maxindexes]
            deviance = np.min(deviances)
            if(deviance <= positive_trigger_length):
                policy_detected_bps += 1
            else:
                file_deviation.append(deviance)
                deviations.append(deviance)

        # Print result
        if verbatim:
            print(f"breakpoint_{i}.txt ->\nExpected breakpoints =\n{breakpoints}")
            print(f"indexes detected :\n{np.array(maxindexes)}")
            if len(file_deviation) > 0 : print(f"added deviance of {file_deviation}")
                
    # Compute median of deviations
    sorted_deviations = np.sort(deviations)
    median_deviation = sorted_deviations[len(sorted_deviations) // 2]

    print("========= RESULTS ========")
    print(f"data directory: {datadir}")
    print(f"Number of matrices tested: {n_dir}")
    print(f"===== Policy ======")
    print(f"True positives = {policy_detected_bps}/{total_bps} ({round(policy_detected_bps/total_bps,2)})")
    print(f"Median deviance for misses = {median_deviation}")


if __name__ == "__main__":
    test_bam_policy()