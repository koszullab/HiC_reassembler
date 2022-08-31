from doctest import testfile
import sys
from pathlib import Path
from tracemalloc import start
import hiscram.reassemble.updating.inversion as upd_inv
from os import mkdir
from Bio import SeqIO, Seq
from Bio.Seq import MutableSeq
import hicstuff.pipeline as hpi
import tensorflow as tf
import keras
from keras.utils.vis_utils import plot_model
from hiscram.genome import Genome
import hicstuff.pipeline as hpi
import numpy as np
from cooler import Cooler
from os import path
import shutil
import pyfastx

fasta = "tutorial_data/W303.fa"
chrom = "chrII"
genome = Genome(pyfastx.Fasta(str(fasta)))
genome.invert(chrom, 124589, 389003)

mod_genome = "test_invert_at_restriction.fa"
with open(mod_genome, "w") as fasta_out:
    seqs = genome.get_seq()
    for chrom, seq in seqs.items():
        fasta_out.write(f">{chrom}\n")
        for frag in seq:
            fasta_out.write(frag)
        fasta_out.write("\n")


hpi.full_pipeline(
    mod_genome,
    "~/Documents/W303_scrambling/reads_chrII_R1.fq.gz",
    "~/Documents/W303_scrambling/reads_chrII_R2.fq.gz",
    aligner="bowtie2",
    out_dir="out_test",
    prefix="test_invert_at_restriction",
    threads=8,
    enzyme=1000,
    mat_fmt="cool",
    force=True,
    no_cleanup=True,
)