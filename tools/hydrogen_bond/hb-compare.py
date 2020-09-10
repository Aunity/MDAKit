#!/usr/bin/env python
import os
import sys
import argparse

import numpy as np
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Compare two msm ensemble hydrogen bond. h1 - h2")
    parser.add_argument("-h1", dest="h1f", help="The msm ensemble hydrogen bond file.",required=True)
    parser.add_argument("-h2", dest="h2f", help="As the h1",required=True)
    parser.add_argument("-o", dest="outf", help="The result file name.", default="hb-msm-comp.txt")

    args = parser.parse_args()

    return args.h1f, args.h2f, args.outf

def main():
    h1f, h2f, outf = parse_args()

    hbd1 = pd.read_csv(h1f, sep="\s+")
    hbd2 = pd.read_csv(h2f, sep="\s+")

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
            col_freq += hbd1_dict[k]["frequence"]
            sfreq.append(hbd1_dict[k]["frequence"])
        else:
            sfreq.append(0)
        if k in hbd2_dict:
            col_atoml = hbd2_dict[k]["label"]
            col_freq -= hbd2_dict[k]["frequence"]
            rfreq.append(hbd2_dict[k]["frequence"])
        else:
            rfreq.append(0)
        freq.append(col_freq)
        atoml.append(col_atoml)

    hb_comp = pd.DataFrame({"label":atoml, "frequence":freq, "key":list(hb_keys),"ref-freq":rfreq,"state-freq":sfreq})


    hb_comp.to_csv(outf, sep=" ", index=None)

if __name__ == "__main__":
    main()
