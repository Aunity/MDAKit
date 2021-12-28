#!/software/anaconda2/envs/msmb/bin/python

import os
import sys
import argparse

import mdtraj as md
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.stats.stats import pearsonr

font = {
        'font.family': 'calibri',
        'font.weight': 'normal',
        'font.size': 20,
        'lines.linewidth':2,
        'figure.dpi':600,
        'figure.autolayout':True
         }
mpl.rcParams.update(font)

def random_coil(res,msm,exp,t,coil_type='v1'):
    '''
    Ref: Wishart DS, Bigam CG, Holm A, Hodges RS, Sykes BD. 1H, 13C and 15N random coil NMR chemical shifts of the common amino acids. I. Investigations of nearest-neighbor effects. J Biomol NMR. 1995 Jan;5(1):67-81. doi: 10.1007/BF00227471. PMID: 7881273.
    '''
    # table 2: RANDOM COIL 13C CHEMICAL SHIFTS FOR THE 20 COMMON AMINO ACIDS WHEN FOLLOWED BY ALANINE 
    random_coilv1 = {
    'CA':{
        "A":52.5, "R":56.0, "N":53.1, "D":54.2,
        "C":58.2, "Q":55.7, "E":56.6, "G":45.1,
        "H":55.0, "I":61.1, "L":55.1, "K":56.2,
        "M":55.4, "F":57.7, "P":63.3, "S":58.3,
        "T":61.8, "W":57.5, "Y":57.9, "V":62.2,
        },
    'CB':{
          "A":19.1, "R":30.9, "N":38.9, "D":41.1,
          "C":28.0, "Q":29.4, "E":29.9, "G":0.0,
          "H":29.0, "I":38.8, "L":42.4, "K":33.1,
          "M":32.9, "F":39.6, "P":32.1, "S":63.8,
          "T":69.8, "W":29.6, "Y":38.8, "V":32.9,
        }
            }
    # table 4: RANDOM COIL 13C CHEMICAL SHIFTS FOR THE 20 COMMON AMINO ACIDS WHEN FOLLOWED BY PROLINE
    random_coilv2 = {
    'CA':{
        "A":50.5, "R":54.0, "N":51.3, "D":52.2,
        "C":56.4, "Q":53.7, "E":54.2, "G":44.5,
        "H":53.3, "I":58.7, "L":53.1, "K":54.2,
        "M":53.3, "F":55.6, "P":61.5, "S":56.4,
        "T":59.8, "W":55.7, "Y":55.8, "V":59.8,
        },
    'CB':{
          "A":18.1, "R":30.2, "N":38.7, "D":40.9,
          "C":27.1, "Q":28.8, "E":29.2, "G":0.0,
          "H":29.0, "I":38.7, "L":41.7, "K":32.6,
          "M":32.4, "F":39.1, "P":30.9, "S":63.3,
          "T":69.8, "W":28.9, "Y":38.3, "V":32.6,
        }
            }
    if coil_type == "v1":
        random_coil = random_coilv1[t]
    else:
        random_coil = random_coilv2[t]

    for i,r in enumerate(res):
        coilValue = random_coil[r]
        exp[i]   -= coilValue
        msm[i]   -= coilValue
    return exp,msm

def parse_arg():
    parser = argparse.ArgumentParser(description='Calculate chemical shift for xtc')
    parser.add_argument('-m', dest='msm', help="msm chemicalshift data",required=True)
    parser.add_argument('-e', dest='exp_f', help="experiment chemicalshift data", default='/home/ymh/mybin/BMR18601.txt')
    parser.add_argument('-o', dest='oup', help="floder to save png, default ='.'", default='.')
    parser.add_argument('-coil_type', help="random_coil data to use, v1 or v2, default:v1", default="v1", choices=['v1','v2'])
    parser.add_argument('-c', dest='Coil', help='Swich. Cal Random coil rather than standard cs.default=False',action="store_false" )
    args = parser.parse_args()
    return args.msm, args.exp_f, args.oup,args.Coil,args.coil_type

def main():
    msm_f,exp_f,rs_p,Coil,coil_type = parse_arg()
    print exp_f
    names = ['CA','CB','N','H']
    data  = pd.read_csv(msm_f, sep='\s+', na_values="****")#, skiprows=1)
    data1 = pd.read_csv(exp_f, sep="\s+", na_values="****")
    for name in names:
        if len(data)!=len(data1):
            print("Two file contains difference rows!!!,exp:msm %d,%d"%(len(data1),len(data)))
            sys.exit()
        fig,ax = plt.subplots(dpi=100, figsize=(6,5))
        comp = pd.DataFrame({'Num':data['Num'], 'Res':data['RES'],'msm':data[name],'exp':data1[name]})
        comp = comp.dropna()
        msm = np.array(comp['msm'])
        exp = np.array(comp['exp'])
        res = np.array(comp['Res'])
        num = comp['Num']

        if Coil and name == "CA":
            exp, msm = random_coil(res,msm,exp,name, coil_type)
            name = 'CA-coil%s'%coil_type

        if Coil and name == "CB":
            exp, msm = random_coil(res,msm,exp,name, coil_type)
            name = 'CB-coil%s'%coil_type
        R,Pv = pearsonr(msm,exp)

        ax_max = np.max(msm)
        ax_min = np.min(msm)
        ax_range = (ax_min,ax_max)
        y = np.arange(ax_min,ax_max+10)
        #fig.set_size_inches(7,6)
        ax.scatter(exp,msm,c='blue')
        ax.plot(y,y,c='red')
        ax.set_xlabel('Expt.(ppm)')
        ax.set_ylabel('Calc.(ppm)')
        ax.set_xlim(ax_range)
        ax.set_ylim(ax_range)
        ax.set_title("%s:%.3f"%(name,R))
        fig.tight_layout()
        fig.savefig(os.path.join(rs_p,'%s.png'%name))

        np.savetxt(os.path.join(rs_p,"%s.txt"%name), np.array([num, res, exp, msm]).T, fmt="%d %s %.4f %.4f")
        plt.close()

if __name__ == "__main__":
    main()
