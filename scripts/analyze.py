import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm

def do_plot(x, dat, fname="fig.png"):
    nlines = dat.shape[1]
    cm_sub = np.linspace(0, 1, nlines)
    colors = [cm.coolwarm(c) for c in cm_sub]

    plt.figure()
    for col_idx in range(nlines):
        plt.plot(x, [1.00] + dat[:,col_idx].tolist() + [0.15], 'o-', color=colors[col_idx])

    plt.xlabel(r"$x'$")
    plt.ylabel(r"$T'$")
    plt.savefig(fname)

def main():
    dat = np.loadtxt("data.txt")
    x = np.linspace(0, 1, 12)
    do_plot(x, dat)

if __name__ == "__main__":
    main()
