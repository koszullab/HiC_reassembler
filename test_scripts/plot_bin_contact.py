import numpy as np
import cooler
import click
import pathlib
from matplotlib import pyplot as plt
from os import path

@click.command()
@click.option(
    "--binsize",
    "-b",
    type=int,
    default=1000
)
@click.argument("dir", type=click.Path(exists=True,path_type=pathlib.Path))
def plot_bin_contact(dir,binsize):
    """
    """
    n = 0
    while path.isdir(dir / "scrambling" / f"scrambled_{n}"):
        n += 1
    
    while True:
        i = int(np.random.uniform(low=0, high=n))
        cool = cooler.Cooler(str(dir / "scrambling" / f"scrambled_{i}" / f"scrambled_{i}.cool"))
        cooler.balance_cooler(cool,store=True)
        mat = cool.matrix(sparse=True, balance=True)

        breakpoints = []
        with open(dir / "scrambling" / f"scrambled_{i}" / f"breakpoint_{i}.txt") as bp:
            for line in bp:
                L = line.strip().split()
                line_breakpoints = L[1].split(':')[1].split('-')
                breakpoints += [int(int(b)/binsize) for b in line_breakpoints]
        breakpoints = np.unique(breakpoints)
        print(breakpoints)
        for bp in breakpoints:
            row = mat[bp].toarray()[0]
            X = np.arange(len(row))
            plt.plot(X,row)
            plt.title(f"{bp}")
            plt.show()
        


if __name__ == "__main__":
    plot_bin_contact()