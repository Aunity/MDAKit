import os
import argparse
import numpy as np
import pandas as pd

import tempfile
from MDAnalysis import Universe
from MDAnalysis.analysis import align


def align_traj(trajfile, topfile, reffile, selstr, outfile):
    if outfile is None:
        outfile = os.path.split(trajfile)[-1][:-4] + '_fit.xtc'
    traj = Universe(topfile, trajfile)
    if os.path.isfile(reffile):
        ref = Universe(reffile)
    else:
        traj.trajectory[reffile]
        tmppdb = tempfile.mktemp(suffix=".pdb")
        grp = traj.select_atoms('all')
        grp.write(tmppdb)
        ref = Universe(tmppdb)
    aligned = align.AlignTraj(traj, ref, filename=outfile, prefix=None, select=selstr)
    aligned.run()
    rmsdfile = outfile[:-4] + '.rmsd'
    np.savetxt(rmsdfile, aligned.rmsd, fmt='%.4f')

def main():
    trajfile = '../example/msm/msm_sample100/set_000000.xtc'
    topfile = '../example/msm/top.pdb'
    align_traj(trajfile, topfile, 10, 'backbone', 'test.xtc')

if __name__ == "__main__":
    main()