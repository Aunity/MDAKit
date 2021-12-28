#!/usr/bin/env python
'''
@Author ymh
@Email  maohuay@hotmail.com
@Date   2021-04-20 14:31:33
@Web    https://github.com/Aunity
'''

import os,sys
import warnings
import tempfile
import numpy as np
import mdtraj as md
import pandas as pd
from multiprocess import Pool
from glob import glob
warnings.filterwarnings('ignore')

# import SESCA modules
SESCA_dir = '/home/ymh/software/SESCA/'
Script_dir = os.path.join(SESCA_dir,"scripts")
sys.path.append(Script_dir)

import SESCA_main as Main
import SESCA_pred as Pred

def predict_CD_spectra_pdb(File, outf):
    SESCA_args = " @pdb %s @write %s @verb 0 @range 150.0,250.0" % (File, outf)
    Processed_Args = Main.Read_Args(SESCA_args.split())
    Data = Main.SESCA_Main(Processed_Args)
    #   Data array 0: SS composition, CD spectra
    #   Data array 1: comparison results (now empty)
    return Data

def predict_CD_spectra_xtc(xtcFile, top, stride=1, outp=None):
    DC = 0
    xtc = md.load(xtcFile, top=top, stride=stride)
    if outp is not None:
        if not os.path.exists(outp):
            os.mkdir(outp)
    else:
        outp = tempfile.mkdtemp()
        DC = 1
    cwd = os.getcwd()
    os.chdir(outp)
    CD_Spectra,SS_Comp = [],[]
    for i in range(len(xtc)):
        pdbFile = 'c%d.pdb' % i
        CDFile = 'c%d.CD' % i
        xtc[i].save_pdb(pdbFile)
        Data = predict_CD_spectra_pdb(pdbFile, CDFile)
        SS_Comp.append(Data[0][0])
        CD_Spectra.append(Data[0][1])
    CD = np.mean(np.array(CD_Spectra), axis=0)
    os.chdir(cwd)
    if DC == 1:
        cmd = 'rm -rf %s' % outp
        os.system(cmd)
    return CD

def func(args):
    xtcFile, top = args
    return predict_CD_spectra_xtc(xtcFile, top)

def obtain_CD_spectra_from_MSM(msmf, trajPath, top, stride=1):
    M = pd.read_pickle(msmf)
    trajfiles = sorted(glob(os.path.join(trajPath, '*.xtc')))
    Wave, CDs, CD = [], [], None
    for trajfile, pi in zip(trajfiles, M.populations_):
        CD_ = predict_CD_spectra_xtc(trajfile, top, stride=stride)
        Wave.append(CD_[:,0])
        CDs.append(CD_[:,1])
        if CD is None:
            CD = CD_[:,1] * pi
        else:
            CD += CD_[:,1] * pi
    return np.c_[Wave[0], CD]

def CD_traj():
    if len(sys.argv[1:]) != 3:
        print('Usage:python %s <traj_file> <topology_file> <out_file>'%sys.argv[0])
        sys.exit(0)
    trajfile, topf, outf = sys.argv[1:]
    data = predict_CD_spectra_xtc(trajfile, topf)
    np.savetxt(outf, data, fmt='%.4f', header='Wave, CD spectra(1000 degrees*cm^2/dmol)')

def CD_pdb():
    if len(sys.argv[1:]) != 2:
        print('Usage:python %s <pdb_file> <out_file>'%sys.argv[0])
        sys.exit(0)
    pdbfile, outf = sys.argv[1:]
    data = predict_CD_spectra_xtc(pdbfile, pdbfile)
    np.savetxt(outf, data, fmt='%.4f', header='Wave, CD spectra(1000 degrees*cm^2/dmol)')

def CD_MSM():
    if len(sys.argv[1:]) != 4:
        print('Usage:python %s <Mobjec_pickle> <sample_folder> <topology_file> <out_file>'%sys.argv[0])
        sys.exit(0)
    msmf, trajPath, top, outf = sys.argv[1:]
    data = obtain_CD_spectra_from_MSM(msmf, trajPath, top, 1)
    np.savetxt(outf, data, fmt='%.4f', header='Wave, CD spectra(1000 degrees*cm^2/dmol)')

if __name__ == '__main__':
    proname = os.path.split(sys.argv[0])[-1]
    if proname == 'CD-pdb':
        CD_pdb()
    elif proname == 'CD-traj':
        CD_traj()
    elif proname == 'CD-msm':
        CD_MSM()
