from random import randint
import numpy as np
from os import path
import click
from pathlib import Path

@click.command()
@click.argument("datadir", type=click.Path(exists=True, path_type=Path))
@click.argument("name", type=str, default="scrambled")
def assemble_dataset_matrix(datadir, name):
    """
    Tool script to assemble a dataset.
    As input, requires a folder containing a dataset.
    The hierarchy of the required dataset corresponds to the one create by the provided scripts. (create_dataset_matrix.py)
    Outputs the imgs.npy and labels.npy required for the matrixdetector.
    """
    scrambledir = datadir / "scrambling"
    normaldir = datadir / "normal"

    imgs = []
    labels = []

    index = 0
    while path.isdir(scrambledir / f"{name}_{index}"):
        try:
            scrambled_imgs = np.load(scrambledir / f"{name}_{index}" / f"{name}_{index}.npy")
            next_scrambled_imgs = np.load(scrambledir / f"{name}_{index}" /  f"next_{name}_{index}.npy")
            for scrambled_img in scrambled_imgs:
                imgs.append(scrambled_img)
                labels.append(1)
            r = randint(0,next_scrambled_imgs.shape[0]-1)
            next_img = next_scrambled_imgs[r]
            imgs.append(next_img)
            labels.append(0)
        except:
            pass
        index += 1

    normal_imgs = np.load(normaldir / "normal_imgs.npy")
    for normal_img in normal_imgs:
        imgs.append(normal_img)
        labels.append(0)

    (datadir / f"dataset").mkdir(parents=True, exist_ok = True)
    imgs = np.array(imgs)
    labels = np.array(labels)
    print(f"imgs shape = {imgs.shape}\nlabels shape = {labels.shape}")
    np.save(datadir / "dataset" / "imgs", imgs)
    np.save(datadir / "dataset" / "labels", labels)

if __name__ == "__main__":
    assemble_dataset_matrix()