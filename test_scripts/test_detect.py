from random import Random, random
from pyparsing import col
from hiscram.detector.matrixdetector import Matrixdetector
import numpy as np
import click
from os import path
import pathlib
import tensorflow as tf

@click.command()
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
@click.argument("datadir", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("name", type=str, default="scrambled")
@tf.autograph.experimental.do_not_convert
def test_detect(datadir,name, chrom,binsize):
    """
    Tests the breakpoints detection on the contact matrix, on a whole dataset. Is design to work on a dataset created by the provided code.
    """

    MatDetect = Matrixdetector()
    MatDetect.load()
    scrambledir = datadir / f"scrambling"

    # Retrieve the number of scrambling matrices available
    n_dir = 0
    while path.isdir(scrambledir / f"{name}_{n_dir}"):
        n_dir += 1

    # Select and detect {number} matrices amongst the dataset
    total_bps = 0
    detected_bps = 0
    false_positives = 0
    policy_detected_bps = 0
    policy_false_positive = 0
    for i in range(n_dir):
        # Select the matrix with number i
        print(i)
        breakpoints = []
        cool_mat = str(scrambledir/ f"{name}_{i}" / f"{name}_{i}.cool")
        # Retreive the breakpoints
        with open(scrambledir / f"{name}_{i}" / f"breakpoints_{i}.txt") as bp:
            for line in bp:
                breakpoints.append(int(int(line)/binsize))

        breakpoints = np.unique(breakpoints)
        total_bps += len(breakpoints)
        # Predict the breakpoints
        inds_SV_detected, probs_SV_detected =  MatDetect.predict(cool_mat, chrom)
        print(f"breakpoint_{i}.txt -> Expected breakpoints = {breakpoints}")
        print(f"indexes detected : {inds_SV_detected}")
        print(f"with probab      : {probs_SV_detected}")
        
        # Measures the true positives VS false negatives
        for bp in breakpoints:
            if bp in inds_SV_detected:
                detected_bps += 1
        
        # Measures the false positives
        for detected in inds_SV_detected:
            if detected not in breakpoints:
                false_positives += 1
        
        # Measures the efficiency of the policy "Take most probable by cluster"
        # Create clusters
        cluster_inds = []
        cluster_probs = []
        for p in range(len(inds_SV_detected)):
            if len(cluster_inds) != 0 and inds_SV_detected[p] == cluster_inds[-1][-1] + 1:
                cluster_inds[-1].append(inds_SV_detected[p])
                cluster_probs[-1].append(probs_SV_detected[p])
            else:
                cluster_inds.append([inds_SV_detected[p]])
                cluster_probs.append([probs_SV_detected[p]])
        # Find most probable by cluster
        most_probables_inds = []
        for c in range(len(cluster_inds)):
            amax = np.argmax(cluster_probs[c])
            most_probable_ind = cluster_inds[c][amax]
            most_probables_inds.append(most_probable_ind)
        print("most probable indexes =", most_probables_inds)
        # Check if the most probable match with the real bp index
        # TP VS FN
        for bp in breakpoints:
            if bp in most_probables_inds: 
                policy_detected_bps += 1
        # FP
        for detected in most_probables_inds:
            if detected not in breakpoints:
                policy_false_positive += 1
    # Print results of analysis
    print("========= RESULTS ========")
    print(f"data directory: {datadir}")
    print(f"Number of matrices tested: {n_dir}")
    print("========== RAW =====")
    print(f"True positives = {detected_bps}/{total_bps} ({round(detected_bps/total_bps,2)})")
    print(f"False positives = {round(false_positives/n_dir,2)} per matrix")
    print("========== Policy ==")
    print(f"True positives = {policy_detected_bps}/{total_bps} ({round(policy_detected_bps/total_bps,2)})")
    print(f"False positives = {round(policy_false_positive/n_dir,2)} per matrix")



if __name__ == "__main__":
    test_detect()