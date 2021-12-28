#!/usr/bin/env python
import os
import sys

import numpy as np
import pandas as pd

from pakg import set_plt

def main():
    print("Usage: python rmsd-plot.py rmsd_lipid.pickl 1.0 rmsd-bilayer.txt")
    rmsdf, cutoff, rmsd_calf = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    rmsd = np.array(pd.read_pickle(rmsdf))
    rmsd_cal = np.loadtxt(rmsd_calf)
    X = np.arange(0, rmsd.shape[1]) * 0.8
    Num = []
    for r in rmsd.T:
        tmp = np.where(r>cutoff)[0]
        Num.append(len(tmp))
    Num = np.array(Num)

    plt = set_plt()
    fig, ax = plt.subplots()
    ax.plot(X, Num, label="Count")
    ax.set_ylim((0, 136))
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("Num.")

    ax2 = ax.twinx()
    ax2.plot(X, rmsd_cal[:,0], label="rmsd", color="red")
    ax2.set_ylim((0.7,1.6))
    ax2.set_ylabel("rmsd (nm)")

    ax2.hlines(cutoff, X[0], X[-1], color="green")
    fig.tight_layout()
    fig.savefig("RMSD.png")

if __name__ == "__main__":
    main()
