#!/usr/bin/env python
import os
import sys
import argparse

import numpy as np
import pandas as pd

def compare_hydrogenbonds(hbfile1:str, hbfile2:str, outfile:str):

    hbd1 = pd.read_csv(hbfile1, sep="\s+")
    hbd2 = pd.read_csv(hbfile2, sep="\s+")

    # hbd1_array = np.array(hbd1)
    # hbd2_array = np.array(hbd2)

    # hbd1_dict = dict(zip(hbd1_array[:,1], hbd1_array))
    # hbd2_dict = dict(zip(hbd2_array[:,1], hbd2_array))

    hbd1_dict = dict([ (irow["key"],irow) for _,irow in hbd1.iterrows()])
    hbd2_dict = dict([ (irow["key"],irow) for _,irow in hbd2.iterrows()])

    hb_keys = set(list(hbd1_dict.keys()) + list(hbd2_dict.keys()))
    freq = []
    atoml= []
    sfreq = []
    rfreq = []
    for k in hb_keys:
        col_atoml = None
        col_freq  = 0
        if k in hbd1_dict:
            col_atoml = hbd1_dict[k]["label"]
            col_freq += hbd1_dict[k]["probability"]
            sfreq.append(hbd1_dict[k]["probability"])
        else:
            sfreq.append(0)
        if k in hbd2_dict:
            col_atoml = hbd2_dict[k]["label"]
            col_freq -= hbd2_dict[k]["probability"]
            rfreq.append(hbd2_dict[k]["probability"])
        else:
            rfreq.append(0)
        freq.append(col_freq)
        atoml.append(col_atoml)

    hb_comp = pd.DataFrame({"label":atoml, "probability":freq, "key":list(hb_keys),"h2-prob":rfreq,"h1-prob":sfreq})


    hb_comp.to_csv(outfile, sep=" ", index=None)

def main():
    parser = argparse.ArgumentParser(description="Compare two msm ensemble hydrogen bond. h1 - h2")
    parser.add_argument("-h1", dest="h1f", help="The ensemble hydrogen bond file.",required=True)
    parser.add_argument("-h2", dest="h2f", help="As the h1",required=True)
    parser.add_argument("-o", dest="outf", help="The result file name.", default="hb-msm-comp.txt")
    args = parser.parse_args()

    compare_hydrogenbonds(args.h1f, args.h2f, args.outf)

if __name__ == "__main__":
    main()