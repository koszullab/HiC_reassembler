#!/usr/bin/env python3
# Script used to generate training data for breakpoint detection on the contact matrix
from argparse import ArgumentError
import pathlib
import click
from os import path

from numpy import require
from create_scrambled_matrix import create_dataset_matrix_scramble
from create_normal_matrix import create_dataset_matrix_normal
from create_full_dataset import create_full_dataset

@click.command()
@click.option(
    "--reads1",
    "-1",
    default=None,
    help="Forward Hi-C reads",
    required=True,
    type=click.Path(exists=True, path_type=pathlib.Path),
)
@click.option(
    "--reads2",
    "-2",
    default=None,
    help="Reverse Hi-C reads",
    required=True,
    type=click.Path(exists=True, path_type=pathlib.Path),
)
@click.option(
    "--binsize",
    "-b",
    default=1000,
    show_default=True,
    help="Size of the bins for the matrices to generate, in base pairs",
)
@click.option(
    "--tmpdir",
    "-t",
    default="./tmp",
    help="Temporary directory to use for the runs",
    type=click.Path(exists=False, path_type=pathlib.Path),
)
@click.option(
    "--name",
    "-n",
    default="scrambled",
    show_default=True,
    type=str,
    help="Name of the scrambled matrices",
)
@click.option(
    "--repeat",
    "-r",
    default=50,
    show_default=True,
    type=int,
    help="Number of matrices to generate",
)
@click.option(
    "--chrom",
    "-c",
    default=None,
    required=True,
    type=str,
    help="Name of the chromosome to scramble",
)
@click.option(
    "--sv_type",
    "-sv",
    default=None,
    type=str,
    help="Type of SVs you want to generate (TRA, INV, DEL). Random between TRA and INV by default."
)
@click.option(
    "--reset",
    "-res",
    is_flag = True,
    help = "Cleans the output directory before running the scrambling."
)
@click.option(
    "--mat_type",
    "-m",
    type=str,
    default="both",
    show_default=True,
    help="Select if you want to generate \"scrambled\" images, \"normal\" images, or \"both\". \"both\" by default will generate the whole dataset. \"scrambled\" will only run some scrambling without saving the dataset (usefull if you want to generate some matrices without regenerating a whole dataset)."
)
@click.option(
    "--aligner",
    "-a",
    type=str,
    show_default=True,
    default="bowtie2",
)
@click.option(
    "--sv_per_matrix",
    "-N",
    type=int,
    default=5,
    show_default=True,
    help="Number of SV to generate per scrambled matrix."
)
@click.option(
    "--wsize",
    "-w",
    type=int,
    show_default=True,
    default=128,
    help="Size of the window to use for image generation."
)
@click.option(
    "--genome_only",
    is_flag = True,
    default=False,
    help="Select if you want to generate only the scrambled genomes without the matrices. (without running hicstuff pipeline)"
)

@click.argument("fasta", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument(
    "outdir", type=click.Path(exists=False, path_type=pathlib.Path)
)
def create_dataset_matrix(
    fasta, outdir, reads1, reads2, binsize, tmpdir, name, repeat, chrom, sv_type, reset, mat_type,aligner,sv_per_matrix,genome_only, wsize
):  
    """
    Creates scrambled matrix aligning reference reads on a scrambled genome.
    By selecting "both" as mat_type, a whole dataset is created and save as npy files for the specified contig.
    By selecting "scrambled", several scrambled matrix can be created with several arguments, without saving any dataset.
    """
    if reads1 is None or reads2 is None:
        raise NotImplementedError(
            "Reads generation not implemented yet. Please provide input reads."
        )
    if chrom is None:
        raise NotImplementedError(
            "Please provide the name of a chromosome to scramble using -c/--chrom"
        )
    if sv_type is not None and sv_type not in ["DEL", "INV", "TRA"]:
        raise ArgumentError(
            "Given sv_type not corresponding to a structural variation. DEL, INS, TRA only are supported."
        )
    if(mat_type not in ["scrambled", "normal", "both"]):
        raise ArgumentError(
            "Wrong mat_type. Only \"scrambled\", \"normal\" or \"both\" allowed."
        )

    if(mat_type == "scrambled"):
        create_dataset_matrix_scramble(fasta, outdir, reads1, reads2, binsize, tmpdir, name, repeat, chrom, sv_type, reset, aligner,sv_per_matrix, genome_only, wsize)
    if(mat_type == "normal"):
        create_dataset_matrix_normal(fasta,outdir,reads1,reads2,binsize,tmpdir,name,repeat,chrom,reset,aligner,genome_only, wsize)
    if(mat_type == "both"):
        create_full_dataset(fasta, outdir, reads1, reads2, binsize, tmpdir, name, repeat, chrom, sv_type, reset, aligner,sv_per_matrix, genome_only, wsize)

if __name__ == "__main__":
    create_dataset_matrix()


        