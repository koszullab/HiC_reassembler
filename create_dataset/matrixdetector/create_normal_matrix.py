#!/usr/bin/env python3
# Script used to generate training data for breakpoint detection on the contact matrix
from argparse import ArgumentError

from sklearn.multiclass import OutputCodeClassifier
from hiscram.scramble import sv
from hiscram.genome import Genome
import hicstuff.pipeline as hpi
import pathlib
import click
import numpy as np
from cooler import Cooler
from os import path
import pyfastx


def create_dataset_matrix_normal(
    fasta, outdir, reads1, reads2, binsize, tmpdir, name, repeat, chrom, reset, aligner, wsize = 128
):  
    """
    Generates images of size wsize extracted of the contact matrix of the given chromosome.
    These images are meant to be given to the CNN as images corresponding to
    positions without a breakpoint.
    The resulting images are saved in "$(outdir)/normal/normal_imgs.npy"
    """

    normaldir = outdir / f"normal"
    hpi.full_pipeline(
        fasta,
        reads1,
        reads2,
        aligner=aligner,
        tmp_dir=str(normaldir / "tmp"),
        out_dir=str(normaldir),
        prefix="normal",
        threads=8,
        enzyme=binsize,
        mat_fmt="cool",
        force=True,
        no_cleanup=True,
    )
    c = Cooler(str(normaldir / f"normal.cool"))
    mat = np.array(c.matrix(balance=False).fetch(chrom))
    N = mat.shape[0]
    h = wsize // 2
    normal_imgs = []
    mat_coords = np.random.randint(low=h, high=N-h,size=repeat)
    for coord in mat_coords:
        win = np.array(mat[(coord - h) : (coord + h), (coord - h) : (coord + h)])
        normal_imgs.append(win)
    normal_imgs = np.array(normal_imgs)
    print(normal_imgs.shape)
    np.save(outdir / f"normal" / "normal_imgs", normal_imgs)