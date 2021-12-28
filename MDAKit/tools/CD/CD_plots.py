#!/usr/bin/env python
'''
@Author ymh
@Email  maohuay@hotmail.com
@Date   2021-04-20 18:36:48
@Web    https://github.com/Aunity
'''

import os,sys
import warnings
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.style.use('~/mybin/ymh.mplstyle')
warnings.filterwarnings('ignore')

def plot_CD(exp, md, png):
    fig, ax = plt.subplots()
    exp = np.array(sorted(exp, key=lambda x: x[0]))
    ax.plot(*exp.T, '--', color='black', linewidth=3)
    ax.plot(md[:,0]-2, md[:,1], color='red', linewidth=3)
    ax.set_xlabel(r'$\mathregular{\lambda}$')
    ax.set_ylabel(r'$\mathregular{[\theta](10^{3}deg*cm^{2}/dmol)}$')
    ax.set_xlim(190,250)
    fig.tight_layout()
    fig.savefig(png, dpi=100)

def main():
    if len(sys.argv[1:]) != 3:
        print('Usage:python %s <exp_CD_spectra> <pre_CD_spectra> <out_png>'%sys.argv[0])
        sys.exit(0)
    expFile, predFile, pngFile = sys.argv[1:]
    exp = np.loadtxt(expFile)
    md = np.loadtxt(predFile)
    plot_CD(exp, md, pngFile)

if __name__ == '__main__':
    main()
