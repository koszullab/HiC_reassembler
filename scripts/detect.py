# Script used to detect SV breakpoints.
from cmath import exp
from email.policy import default
import pathlib
import numpy as np

from os import mkdir
import click
from hiscram.detector.bamdetectormax import Bamdetectormax
from hiscram.detector.matrixdetector import Matrixdetector

def is_close(position, tab, trigger):
    """
    Tool function returning True if the given position is closer than trigger to
    one of the element of the tab.
    """
    for elem in tab:
        if abs(position-elem) < trigger:
            return True
    return False 

@click.command()
@click.option(
    "--binsize", 
    "-b",
    type=int,
    default=1000,
    help="Binsize used to create the Hi-C matrix.",
)
@click.option(
    "--matrix_only",
    is_flag=True,
    default=False,
    help="In order to stop the detection after the matrix detection."
)
@click.option(
    "--tmpdir",
    "-T",
    type=click.Path(exists=False),
    default="./tmpdir",
)
@click.argument("chrom_name", type=str)
@click.argument("cool", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("for_bam", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("rev_bam", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("sequence", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("breakpoints", type=click.Path(exists=True, path_type=pathlib.Path), required=False, default=None)
def detect(binsize, matrix_only,tmpdir, chrom_name, cool, for_bam, rev_bam, sequence, breakpoints):
    """
    Executes the detection of the breakpoints. First, on the matrix at bin resolution, and then
    at base pair resolution with alignment file.
    The code is suited for the output of hicstuff pipeline (no cleanup).
    """

    # Create temporary drectory
    try:
        mkdir(tmpdir)
    except:
        pass

    # Detection on Hi-C
    print("##### Detection on Matrix at bin level #####")
    MatDetect = Matrixdetector()
    MatDetect.load()
    inds_SV_detected, probs_SV_detected = MatDetect.predict(str(cool), chrom_name)
    print("mat index: ", inds_SV_detected)
    
    if breakpoints is not None:
        expected = []
        with open(breakpoints) as bp:
            for line in bp:
                L = line.strip()
                expected += [int(L)]
        bin_expected = [int(bp/binsize) for bp in expected]
        found = 0
        for bp in inds_SV_detected:
            if bp in bin_expected:
                found += 1
        print("expected    :", bin_expected)
        print(f"accuracy = {found}/{len(bin_expected)}")

    if matrix_only: return

    # Detection on BAM
    print("##### Detection on Alignment file at base pair level #####")
    BamDetect = Bamdetectormax(chrom_name,binsize,str(for_bam), str(rev_bam), ref_seq=str(sequence))
    predictions = BamDetect.predict()
    print("predictions :", np.sort(predictions))
    
    if breakpoints is not None:
        expected = np.sort(np.unique(expected))
        print("expected    :", expected)
        found = 0
        for pos in expected:
            if is_close(pos, predictions, 10):
                found += 1
        print(f"accuracy (+-10) = {found}/{len(expected)} ({int((found/len(expected))*100)}%)")

    # Save detected bps
    with open("data/output/detected_breakpoints.txt", 'w') as detected_bps:
        for bp in predictions:
            detected_bps.write(f"{bp}\n")

if __name__ == "__main__":
    detect()

def is_close(position, tab, trigger):
    """
    Tool function returning True if the given position is closer than trigger to
    one of the element of the tab.
    """
    for elem in tab:
        if abs(position-elem) <= trigger:
            return True
    return False 