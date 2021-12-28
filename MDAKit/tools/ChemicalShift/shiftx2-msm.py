#!/software/anaconda2/envs/msmb/bin/python

import os
import sys
import argparse

import mdtraj as md
import pandas as pd
import numpy as np

#from msmbuilder.io import load_meta, save_meta
from multiprocessing import Pool


def job(cmd):
    os.system(cmd)

def shift_xtc(xtc_p, T=1):
    xtcs = [os.path.join(xtc_p,xtc) for xtc in os.listdir(xtc_p) if xtc.endswith('.xtc')]
    cmds = []
    for xtc in xtcs:
        cmd = "shiftx2-xtc.py -x %s"%(xtc)
        cmds.append(cmd)
    #os.system(cmd)
    p = Pool(T)
    p.map(job,cmds)

def parse_arg():
    parser = argparse.ArgumentParser(description='Calculate chemical shift for xtc')
    parser.add_argument('-c', dest='csp', help="cs result by shiftx2", required=True)
    parser.add_argument('-p', dest='msm', help="msm pikle data M.pikle or pcca.pikle")
    parser.add_argument('-o', dest='oup', help="result file name, default=chemicalshift.txt", default='chemicalshift.txt')
    parser.add_argument('-T', dest='tlp', help="number of threads", default=1, type=int)
    args = parser.parse_args()
    return args.csp, args.msm, args.oup, args.tlp

def main():
    cs_p,msm,rs_p,T = parse_arg()
    M = pd.read_pickle(msm)
    if hasattr(M,'stationary_distribution'):
        P = M.stationary_distribution
        css = [(os.path.join(cs_p,cs),int(cs.split('.')[0][4:])) for cs in sorted(os.listdir(cs_p)) if cs.endswith('.cs')]
    else:
        P = M.populations_
        css = [(os.path.join(cs_p,cs),int(cs.split('.')[0])) for cs in sorted(os.listdir(cs_p)) if cs.endswith('.cs')]
    if hasattr(M,'n_macrostates'):
        macro_p = []
        n = M.n_macrostates
        for i in range(n):
            index = np.where(M.microstate_mapping_==i)
            macro_p.append(P[index])
        total = [sum(m) for m in macro_p]
        P = np.array(total)
    n = len(P)
    if len(css)!=n:
        print(len(P), len(css))
        raise Exception("Error! number of state not equal to the xtc's number")
        sys.exit(1)
    data_msm = None
    index = ['Num','CA','CB','CO','N','H','HA']
    idinfo = None
    for i,(cs,si) in enumerate(css):
        si = i
        data = pd.read_csv(cs, na_values="****", sep=" ")
        if idinfo is None:
            idinfo = data[["Num","RES"]]
        data = data[index]
        #print(data.shape)
        p = [1.0/n]
        p.extend([P[si]]*(data.shape[1]-1))
        data *= p
        if data_msm is None:
            data_msm = data
        else:
            data_msm += data
    data_msm['Num'] = idinfo['Num']
    rsd = pd.merge(idinfo, data_msm, how='right', on=['Num'])
    rsd.to_csv(rs_p, na_rep='****',float_format="%.2f",index=False, sep=" ")
    
if __name__ == "__main__":
    main()
