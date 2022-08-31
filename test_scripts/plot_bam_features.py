from hiscram.detector.bam_functions import (
    bam_region_coverage, bam_region_read_ends, parse_ucsc_region)

import click
from matplotlib import pyplot as plt
import numpy as np
import pathlib
import pysam as ps
import typing

@click.command()
@click.option(
    "--region",
    "-r",
    type=str,
    required=True,
    help="Region to test",
)
@click.option(
    "--mean",
    "-mean",
    is_flag=True,
    help="Activate this flag to mean the signal",
)
@click.argument("bam", type=click.Path(exists=True,path_type=pathlib.Path))
def plot_bam_features(bam,region,mean):
    """
    Tool script that plots several coverages of reads (end+start, clipped) for the given region. 
    """
    chrom, start, end = parse_ucsc_region(region)

    
    start_arr, end_arr = bam_region_read_ends(str(bam), region, side="both", clipped=False)
    if mean:
        start_arr = mean3(start_arr)
        end_arr = mean3(end_arr)
    sum_array = start_arr + end_arr
    start_arr_clip, end_arr_clip = bam_region_read_ends(str(bam), region, side="both", clipped=True)
    if mean:
        start_arr_clip = mean3(start_arr_clip)
        end_arr_clip = mean3(end_arr_clip)
    
    coverage = bam_region_coverage(str(bam),region)
    sum_array_clip = start_arr_clip + end_arr_clip
    X = np.arange(start,end)
    sum_divided = []
    sum_clip_divided = []
    percent_clipped = []
    for i in range(len(X)):
        sum_divided.append(sum_array[i]/ max(coverage[i],1))
        sum_clip_divided.append(sum_array_clip[i]/ max(coverage[i],1))
        percent_clipped.append(sum_array_clip[i]**2 / max(sum_array[i],1))

    fig, ax = plt.subplots(5,1)

    # Sum of start+end of all aligned reads in the region
    ax[0].plot(X,sum_array)
    ax[0].set_ylabel("Start+end")

    # Sum divided by coverage
    ax[1].plot(X,sum_divided)
    ax[1].set_ylabel("(Start+end)/coverage")

    # Sum of start+end of clipped aligned reads in the region
    ax[2].plot(X,sum_array_clip)
    ax[2].set_ylabel("clipped(start+end)")

    # divided by coverage
    ax[3].plot(X,sum_clip_divided)
    ax[3].set_ylabel("clipped(start+end)/coverage")

    # percentage of clipped reads * clipped
    ax[4].plot(X,percent_clipped)
    ax[4].set_ylabel("clipped(start+end)**2/(Start+end)")
    ax[4].set_xlabel(f"max : {np.argmax(percent_clipped) + start}")

    plt.show()

def mean3(arr):
    mean_arr = (
            arr
            + np.concatenate((arr[1:], np.zeros(1)))
            + np.concatenate((np.zeros(1), arr[: len(arr) - 1]))
        ) // 3
    assert(len(mean_arr) == len(arr))
    return mean_arr

if __name__ == "__main__":
    plot_bam_features()