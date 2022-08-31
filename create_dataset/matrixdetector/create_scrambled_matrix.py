#!/usr/bin/env python3
# Script used to generate training data for breakpoint detection on the contact matrix
from hiscram.genome import Genome
import hicstuff.pipeline as hpi
import numpy as np
from cooler import Cooler
from os import path
import shutil
import pyfastx

def create_dataset_matrix_scramble(
    fasta, outdir, reads1, reads2, binsize, tmpdir, name, repeat, chrom, sv_type, reset, aligner,sv_per_matrix, genome_only, wsize
):
    """
    Generate scrambled matrixes in the specified output directory.
    Each scrambling corresponds to 1 folder, in which the result of hicstuff pipeline
    is saved after the given genome is scrambled.
    """

    outdir.mkdir(parents=True, exist_ok=True)
    # To avoid overwriting already existing scramblings. Compute starting index
    # for the new folders.
    starting_index = 0
    if not reset:
        while path.isdir(outdir / f"scrambling/{name}_{starting_index}"):
            starting_index += 1
    else:
        shutil.rmtree(outdir / f"scrambling")
    
    hsize = int(wsize/2)

    for i in range(starting_index, starting_index + repeat):
        # Load genome
        genome = Genome(pyfastx.Fasta(str(fasta)))
        for _ in range(sv_per_matrix):
            index = 0

            # Determine which SV to create. Random btw INV and TRA if not specified by the user.
            if sv_type is None:
                sv_tag = np.random.choice(["INV", "TRA"])
            else:
                sv_tag = sv_type
            
            # Length of the SV to create. These constants are OK, but it can be changed.
            N = len(genome.chroms[chrom])
            sv_size = int(np.random.normal(N/10, N/20))
            sv_size = max(sv_size, 20*binsize)
            print("sv_size = ", sv_size)
            start = np.random.randint(hsize*binsize,N - sv_size - hsize*binsize)
            end = start + sv_size
            
            # Generate the sv
            if sv_tag == "INV":
                genome.invert(chrom, start, end)
            if sv_tag == "DEL":
                genome.delete(chrom, start, end)
            if sv_tag == "TRA":
                insert_pos = start
                while insert_pos >= start and insert_pos <= end:
                    insert_pos = np.random.randint(hsize*binsize,N - sv_size - hsize*binsize)
                # inv = np.random.choice([False, True]) Not yet compatible
                inv = False
                print(f"{start}-{end} -> {insert_pos}")
                # Interchromosomal translocation not yet implemented
                genome.translocate(chrom,insert_pos,chrom,start,end,inv)

        # Delete empty fragments
        genome.chroms[chrom].clean_frags()
        # This is for debugging purpose. I print the resulting fragments of the genome.
        index = 0
        for frag in genome.chroms[chrom].frags:
            print(frag, " at ", index)
            index += len(frag)
        
        # Create the output directory for this scrambling
        scramblingdir = f"scrambling/{name}_{i}"
        (outdir / scramblingdir).mkdir(parents=True, exist_ok=True)
        mod_genome = str(outdir / scramblingdir / f"mod_genome_{i}.fa")
        
        # Save the breakpoints using the boundaries of the fragments. This is the right way.
        breakpoints = genome.chroms[chrom].boundaries[1:-1]
        with open(outdir / scramblingdir / f"breakpoints_{i}.txt", 'w') as bp_file:
            for bp in breakpoints:
                bp_file.write(f"{bp}\n")
            
        # Save the scrambled genome
        with open(mod_genome, "w") as fasta_out:
            seqs = genome.get_seq()
            for chrom, seq in seqs.items():
                fasta_out.write(f">{chrom}\n")
                for frag in seq:
                    fasta_out.write(frag)
                fasta_out.write("\n")

        if genome_only:
            continue
        
        # Run hicstuff pipeline on the scrambled genome
        name_scrambling = f"{name}_{i}"
        hpi.full_pipeline(
            mod_genome,
            reads1,
            reads2,
            aligner=aligner,
            tmp_dir=str(outdir / scramblingdir / "tmp"),
            out_dir=str(outdir / scramblingdir),
            prefix=name_scrambling,
            threads=8,
            enzyme=binsize,
            mat_fmt="cool",
            force=True,
            no_cleanup=True,
        )

        # Retrieve the contact matrix with cooler
        c = Cooler(str(outdir / scramblingdir / f"{name_scrambling}.cool"))
        mat = np.array(c.matrix(balance=False).fetch(chrom))

        # Array containing the coords of the bp on the matrix
        mat_coords = []
        for bp in breakpoints:
            mat_coords.append(int(bp/binsize))
        mat_coords = np.unique(mat_coords)

        sv_imgs = []
        next_sv = []

        # Fetch the image to be given to the CNN.
        # Image of size wsize (128 by default), centered around the breakpoints
        for coord in mat_coords:
            win = np.array(mat[(coord - hsize) : (coord + hsize), (coord - hsize) : (coord + hsize)])
            sv_imgs.append(win)

            # Compute the images right next to the breakpoints.
            # This is supposed to increase the precision of the CNN, but is probably a bad idea.
            rand_offset = np.random.randint(-2,3,size=2)
            for r in rand_offset:
                try:
                    win_next = np.array(mat[(coord - hsize + r) : (coord + hsize + r), (coord - hsize + r) : (coord + hsize + r)])
                except:
                    continue
                next_sv.append(win_next)
        
        # Finally, save the images as .npy files
        np.save(outdir / scramblingdir / f"{name_scrambling}", sv_imgs)
        np.save(outdir / scramblingdir / f"next_{name_scrambling}", next_sv)
        