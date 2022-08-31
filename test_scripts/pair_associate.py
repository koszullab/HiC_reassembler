from turtle import color
import pypairix
import click
import pathlib
import os
from matplotlib import pyplot as plt
import numpy as np

@click.command()
@click.option(
    "--read_length",
    "-rl",
    type=int,
    default=36,
)
@click.option(
    "--winsize",
    "-w",
    type=int,
    default=500,
    help="Length of the query for each side of the given breakpoint",
)
@click.argument("pairfile", type=click.Path(exists=True,path_type=pathlib.Path))
@click.argument("chrom", type=str)
@click.argument("bp", type=int)
def pair_associate(pairfile,bp,read_length,winsize,chrom):
    """
    Plots the positions of the pair reads associated with the reads mapping next to te given position.
    """

    # Gzip the pairfile, because pypairix works only with a gzip file
    os.system(f"bgzip {pairfile} -c > {pairfile}.gz")
    os.system(f"pairix {pairfile}.gz")
    pairs_data = pypairix.open(f"{pairfile}.gz")

    # Sequence where we want the reads to be associated from.
    query_left = f"{chrom}:{bp-winsize}-{bp}"
    query_right = f"{chrom}:{bp}-{bp+winsize}"

    # Queries for the pair file
    it_for_left = pairs_data.querys2D(f"{query_left}|*")
    it_rev_left = pairs_data.querys2D(f"*|{query_left}")
    
    it_for_right = pairs_data.querys2D(f"{query_right}|*")
    it_rev_right = pairs_data.querys2D(f"*|{query_right}")

    # Create coverage of pairs with a dictionary
    # Coverage on the left of the bp
    coverage_left = {}

    for read in it_for_left:
        print(read)
        pair_loc = int(read[4])
        sign = read[6]
        for i in range(read_length):
            if sign == '+':
                index = pair_loc + i
            else:
                index = pair_loc - i
            if coverage_left.get(index) is None:
                coverage_left[index] = 0
            coverage_left[index] += 1

    for read in it_rev_left:
        print(read)
        pair_loc = int(read[2])
        sign = read[5]
        for i in range(read_length):
            if sign == '+':
                index = pair_loc + i
            else:
                index = pair_loc - i
            if coverage_left.get(index) is None:
                coverage_left[index] = 0
            coverage_left[index] += 1
    print(coverage_left)

    # Coverage on right of the bp
    coverage_right = {}
    for read in it_for_right:
        print(read)
        pair_loc = int(read[4])
        sign = read[6]
        for i in range(read_length):
            if sign == '+':
                index = pair_loc + i
            else:
                index = pair_loc - i
            if coverage_right.get(index) is None:
                coverage_right[index] = 0
            coverage_right[index] += 1

    for read in it_rev_right:
        print(read)
        pair_loc = int(read[2])
        sign = read[5]
        for i in range(read_length):
            if sign == '+':
                index = pair_loc + i
            else:
                index = pair_loc - i
            if coverage_right.get(index) is None:
                coverage_right[index] = 0
            coverage_right[index] += 1
    print(coverage_right)

    # Plot the graphs
    for k in coverage_left.keys():
        plt.vlines(k, 0, coverage_left[k], color='b')
    for k in coverage_right.keys():
        plt.vlines(k, 0, -coverage_right[k], color='r')
    plt.show()

if __name__ == "__main__":
    pair_associate()