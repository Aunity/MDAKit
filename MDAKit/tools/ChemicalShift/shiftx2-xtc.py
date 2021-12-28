#!/software/anaconda2/bin/python

import os
import sys
import argparse

import mdtraj as md
import pandas as pd
import numpy as np
import tempfile

def average():
    css = [cs for cs in os.listdir('.') if cs.endswith('.cs')]
    names = ['Num','RES','CA','CB','CO','N','H','HA']
    dtypes = {'Num':int, "RES":str, "CA":np.float64,"CB":np.float64,
              'CO':np.float64,'N':np.float64,'H':np.float64,'HA':np.float64,}
    index = ['Num','CA','CB','CO','N','H','HA']

    Num = None
    Res = None

    n_cs = len(css)
    data_avg = None
    for cs in css:
        data0 = pd.read_table(cs,skiprows=2, sep="\s+", names=names, na_values='****', dtype=dtypes)
        data = data0[index]
        if Num is None:
            Num = np.array(data0['Num'])
        if Res is None:
            Res = np.array(data0['RES'])
        if data_avg is None:
            data_avg = data
        else:
            data_avg += data
    data_avg /= len(css)
    idinfo = pd.DataFrame({"Num":Num,"RES":Res})
    data_avg['Num'] = idinfo['Num']
    rsd = pd.merge(idinfo,data_avg,how='left',on=['Num'])
    return rsd

def parse_arg():
    parser = argparse.ArgumentParser(description='Calculate chemical shift for xtc')
    parser.add_argument('-x', dest='xtc', help="xtc file for cal chemical shift", required=True)
    parser.add_argument('-t', dest='top', help="top file for the xtc file. Default=~/GGBP/2fvy.pdb", default='/home/ymh/GGBP/2fvy.pdb')
    parser.add_argument('-p', dest="ph", help="PH value, default PH=7.0", default=7.0, type=float)
    parser.add_argument('-Tem', dest='Tem', help="Temperature, default Tem=300", default=310, type=int)
    parser.add_argument('-shiftx1', help="Version of shiftx1", action="store_true")

    args = parser.parse_args()
    return args.xtc, args.top, args.ph, args.Tem, args.shiftx1

def main():
    xtc_p,top,PH,Tem,shiftx1 = parse_arg()

    if shiftx1:
        shiftx1 = "-x"
    else:
        shiftx1 = ""
    if xtc_p.endswith('.pdb'):
        xtcs = md.load_pdb(xtc_p)
    else:
        xtcs = md.load_xtc(xtc_p, top)
    pwd = os.getcwd()
    shift_temp = "shift_temp"
    if not os.path.exists(shift_temp):
        os.mkdir(shift_temp)
    rs_p = xtc_p.split('/')[-1]
    temp = rs_p[:-4]
    rs_p = rs_p +'.cs'
    if not os.path.exists(temp):
        os.mkdir(temp)
    os.chdir(temp)

    for i in range(len(xtcs)):
        xtcs[i].save_pdb("%d.pdb"%i)
    tempf = tempfile.mkdtemp()
    cmd = "shiftx2.py -a BACKBONE -b '*.pdb' -f TABULAR -p %f -t %d -z %s %s"%(PH,Tem,tempf,shiftx1)
    os.system(cmd)
    data_avg = average()
    #print(data_avg[:3])
    data_avg.to_csv(rs_p, na_rep='****',float_format="%.2f",index=False, sep=" ")
    os.system("mv %s %s"%(rs_p,os.path.join(pwd,shift_temp)))
    os.chdir(pwd)
    os.system("rm -rf %s"%temp)

if __name__ == "__main__":
    main()

