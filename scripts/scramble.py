#!/usr/bin/env python3
# Script used to generate training data for the model
# Will scramble a genome and generate a Hi-C matrix from it.
from hiscram.scramble import sv
import hicstuff.pipeline as hpi
import pathlib
import click


@click.command()
@click.option(
    "--reads1",
    "-1",
    default=None,
    help="Forward Hi-C reads",
    type=click.Path(exists=True, path_type=pathlib.Path),
)
@click.option(
    "--reads2",
    "-2",
    default=None,
    help="Reverse Hi-C reads",
    type=click.Path(exists=True, path_type=pathlib.Path),
)
@click.option(
    "--binsize",
    "-b",
    default=2000,
    show_default=True,
    help="The resolution of matrices to generate, in basepair",
)
@click.option(
    "--tmpdir",
    "-t",
    default="./tmp",
    help="Temporary directory to use for the runs.",
    type=click.Path(exists=False, path_type=pathlib.Path),
)
@click.option(
    "--profile",
    "-p",
    default="Yeast",
    type=str,
    help="SV size and frequency profile to use. Can be one of: Yeast, Human, Debug",
)
@click.argument("fasta", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument(
    "outdir", type=click.Path(exists=False, path_type=pathlib.Path)
)
def run_scrambles(fasta, outdir, reads1, reads2, binsize, profile, tmpdir):
    """
    This is the orchestrator function that handles the end-to-end pipeline. For
    each scramble run, it will:
    0. Select a random region of a random chromosome
    1. Edit the selected region to add SV
    2. Realign the reads and generate a (cool) Hi-C map from the scrambled genome
    3. Extract windows around each SV as well as random (negative) windows
    4. Store the windows and associated labels into npy files.
    5. Store the whole matrix before and after SV.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    if reads1 is None or reads2 is None:
        raise NotImplementedError(
            "Reads generation not implemented yet. Please provide input reads."
        )

    # Generate random structural variations and apply them to the genome
    mixer = sv.Mixer(str(fasta), config_profile=profile)
    mixer.generate_sv()
    mod_genome = str(outdir / "mod_genome.fa")
    with open(mod_genome, "w") as fa:
        mixer.write_genome(fa)
    with open(outdir / "breakpoints.bedpe", "w") as bp:
        mixer.write_breakpoints(bp)
    # Generate contact map using the edited genome
    hpi.full_pipeline(
        mod_genome,
        reads1,
        reads2,
        aligner="bowtie2",
        tmp_dir=tmpdir,
        out_dir=outdir,
        prefix="scrambled",
        threads=8,
        enzyme=binsize,
        mat_fmt="cool",
    )


if __name__ == "__main__":
    run_scrambles()  # pylint: disable=no-value-for-parameter
