# Script used to reassemble a matrix, given the breakpoints.
import pathlib
import numpy as np
from Bio import SeqIO

from os import mkdir
import click
from hiscram.detector.pair_combiner import PairCombiner
import hicstuff.pipeline as hpi
from hiscram.reassemble.reassembler import Reassembler

@click.command()
@click.option(
    "--binsize",
    "-b",
    type=int,
    default=1000,
    show_default=True,
    help="Binsize used to create the Hi-C matrix.",
)
@click.option(
    "--readlength",
    "-rl",
    type=int,
    show_default=True,
    default=36,
)
@click.option(
    "--tmpdir",
    "-T",
    type=click.Path(exists=False),
    default="./tmpdir",
)
@click.option(
    "--outdir",
    "-o",
    type=click.Path(exists=True),
    default="data/output/",
)
@click.argument("chrom_name", type=str)
@click.argument("cool", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("pairfile", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("ref_seq", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("breakpoints", type=click.Path(exists=True, path_type=pathlib.Path))
def reassemble(binsize,readlength,chrom_name,cool,pairfile,breakpoints,ref_seq,tmpdir,outdir):
    """
    Given a list of knows breakpoints, executes the SV assembly and the reassembly.
    """

        # Create temporary drectory
    try:
        mkdir(tmpdir)
    except:
        pass
    
    print("##### SV assembly using pair file #####")
    SVCombiner = PairCombiner(chrom_name, breakpoints, pairfile, binsize, readlength, tmpdir)
    info_sv = SVCombiner.combine()

    print("##### Genome reassembly #####")
    reassembler = Reassembler(info_sv,cool,ref_seq,chrom_name,binsize)
    mat_reassembled, seq_reassembled = reassembler.reassembly()

    # Write reassembled genome
    try:
        mkdir(outdir)
    except:
        pass
    np.save(outdir + "mat_reassembled.npy", mat_reassembled)
    with open(outdir + "seq_reassembled_final.fa", "w") as fa_out:
        rec = SeqIO.SeqRecord(seq=seq_reassembled, id=chrom_name, description="")
        SeqIO.write(rec, fa_out, format="fasta")


if __name__ == "__main__":
    reassemble()
