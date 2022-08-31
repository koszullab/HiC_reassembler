# Runs the entire pipeline detection + association + reassembly
from os.path import join
import pathlib
import numpy as np
from Bio import SeqIO

from os import mkdir
import click
from hiscram.detector.pair_combiner import PairCombiner
from hiscram.detector.bamdetectormax import Bamdetectormax

from hiscram.detector.matrixdetector import Matrixdetector
from hiscram.detector.pair_combiner import PairCombiner
import hicstuff.pipeline as hpi
from hiscram.reassemble.reassembler import Reassembler

@click.command()
@click.option(
    "--binsize", "-b",type=int, default=1000, help="Binsize used to create the Hi-C matrix.",
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
    type=click.Path(exists=False),
    default="./tmpdir",
)
@click.option(
    "--outdir",
    "-o",
    type=click.Path(),
    default="data/output/",
)
@click.argument("chrom_name", type=str)
@click.argument("cool", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("for_bam", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("rev_bam", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("pairfile", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("ref_seq", type=click.Path(exists=True, path_type=pathlib.Path))
def pipeline(binsize,readlength,chrom_name,cool, for_bam, rev_bam, pairfile, ref_seq, tmpdir, outdir):
    """
    Executes the entire pipeline for the reassembly :
    - Detection on the contact map with the CNN,
    - Detection at base pair precision using the alignement files,
    - Breakpoint association using the pair files
    - Reassembly using the regression model to determine the reassembly order.

    """
    # Create temporary drectory
    try:
        mkdir(tmpdir)
        mkdir(outdir)
    except:
        pass


    # Detection on Hi-C
    print("##### Detection on Matrix at bin level #####")
    MatDetect = Matrixdetector()
    MatDetect.load()
    inds_SV_detected, probs_SV_detected = MatDetect.predict(str(cool), chrom_name)
    print("mat index: ", inds_SV_detected)

    # Detection on BAM
    print("##### Detection on Alignment file at base pair level #####")
    BamDetect = Bamdetectormax(chrom_name,binsize,str(for_bam), str(rev_bam), ref_seq=str(ref_seq))
    predictions = BamDetect.predict()

    print("predictions :", np.sort(predictions))
    
    print("##### SV assembly using pair file #####")
    SVCombiner = PairCombiner(chrom_name, predictions, pairfile, binsize, readlength, tmpdir)
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
    pipeline()