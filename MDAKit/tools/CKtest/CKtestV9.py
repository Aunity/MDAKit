#!/usr/bin/env python
import os
import sys

import tempfile
import argparse

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from multiprocessing import Pool
mpl.style.use("~/mybin/ymh.mplstyle")

def CKtestRound_subjob(args):
    n_states = args[0]
    indexes = args[1]
    lagstep = args[2]
    vartype = args[3]
    _transitionsDICT = args[4]
    counts = np.zeros((n_states, n_states), dtype=float)
    _transitions = []

    for i in indexes:
        _transitions.append(np.row_stack(_transitionsDICT[i]))

    transitions = np.hstack(_transitions)
    #print("debug")
    #print(transitions.shape)
    C = coo_matrix((np.ones(transitions.shape[1], dtype=int), transitions),
                   shape=(n_states, n_states))
    counts = counts + np.asarray(C.todense())

    counts_sum = counts.sum(axis=1)
    counts_sum[counts_sum==0] = 1.0
    counts_i = counts.diagonal()
    TP = counts_i/counts_sum
    return TP

def _transition_counts(sequences, lag_time):
    none_to_nan = np.vectorize(lambda x: np.nan if x is None else x,
                                           otypes=[np.float])
    neg_to_nan = np.vectorize(lambda x: np.nan if x < 0 else x,
                                           otypes=[np.float])
    try:
        classes = np.unique(np.concatenate(sequences))
    except:
        print("sequence should not contains None elements!")
        sys.exit(0)
    #classes = np.array(set(np.concatenate(sequences)))
    contains_nan = (classes.dtype.kind == 'f') and np.any(np.isnan(classes))
    contains_none = any(c is None for c in classes)
    contains_mi = any(c<0 for c in classes)
    _transitionsDICT = {}

    for i, y in enumerate(sequences):
        y = np.asarray(y)
        from_states = y[: -lag_time: 1]
        to_states = y[lag_time::1]

        if contains_none:
            from_states = none_to_nan(from_states)
            to_states = none_to_nan(to_states)

        if contains_mi:
            from_states = neg_to_nan(from_states)
            to_states = neg_to_nan(to_states)
        if contains_nan or contains_none or contains_mi:
            # mask out nan in either from_states or to_states
            mask = ~(np.isnan(from_states) + np.isnan(to_states))
            from_states = from_states[mask]
            to_states = to_states[mask]
        _transitionsDICT[i] = ((from_states, to_states))

    return _transitionsDICT


def CKtestRound(nc, lagstep, labels, Thread=1, random_times=200):
    nt = len(labels);   # number of trajectorys
    random_TPs = np.zeros((nc,random_times))
    var_bar =0.50

    indexes = []
    labels = np.array(labels)

    for i in range(random_times):
        index = np.random.choice(np.arange(nt),int(nt*var_bar))
        indexes.append(index)
    _transitionsDICT = _transition_counts(labels, lagstep)
    #fname = tempfile.mktemp()
    #pd.to_pickle(_transitionsDICT, fname)
    args = [ (nc, index, lagstep, False, _transitionsDICT) for label in indexes ]

    p = Pool(Thread)

    #CKtestRound_subjob(args[0])
    random_TPs = np.array([CKtestRound_subjob(arg)  for arg in args]).T
    #random_TPs = np.array(p.map(CKtestRound_subjob, args)).T

    #TP = random_TPs.mean(axis=1)
    #var = random_TPs.var(axis=1)
    TP  = CKtestRound_subjob((nc, range(nt), lagstep, False, _transitionsDICT))
    print(random_TPs, TP.reshape((nc,1)))
    var = np.sqrt(np.sum((random_TPs-TP.reshape((nc,1)))**2/float(random_times), axis=1))
    del _transitionsDICT
    return TP, var


class ChapmanKolmogorovTest(object):

    def __init__(self, msm, dt, du, outp, nT, nstates, lagtime, pi, nr=3):# lagtime,unit,nc,figp,dt,n_thread,P,nround=3,seqlens=None):
        self.msm = msm
        self.outp = outp
        self.nround = nr
        self.nstates = nstates
        self.lagtime = lagtime
        self.pi = pi
        self.dt = dt
        self.unit = du
        self.MDTP = None
        self.VAR = None
        self.MSMTP = None

        self.n_thread = nT
        self._flag = False

        if not os.path.exists(self.outp):
            os.mkdir(self.outp)

        self.MSMTPf = os.path.join(self.outp, "ProbMSM")
        self.MDTPf  = os.path.join(self.outp, "ProbMD")
        self.VARf    = os.path.join(self.outp, "var")

    def CK_MD(self,labels):
        '''
        MD-left
        '''
        self.labels = labels
        nc = self.nstates
        nround = self.nround
        lagtime = self.lagtime

        if os.path.exists(self.MDTPf) and os.path.exists(self.VARf):
            self.MDTP,self.VAR = np.loadtxt(self.MDTPf, dtype=float), np.loadtxt(self.VARf, dtype=float)
            return

        TPs = np.ones((nc,nround+1))
        Var = np.zeros((nc,nround+1))
        for i in range(1, nround+1):
            print("CK-MD round %d/%d"%(i+1, nround+1))
            lagstep = i * lagtime
            print(i,TPs.shape)
            TPs[:,i],Var[:,i] = CKtestRound(self.nstates, lagstep, labels, self.n_thread)
        self.MDTP,self.VAR = TPs, Var
        self._flag = True

    def CK_MSM(self,T=None):
        '''
        T -transition_mat
        '''
        if T is None:
            T = self.msm.transmat_

        nc, nround = self.nstates, self.nround
        TPs = np.ones((nc,nround+1))

        for j in range(nc):
            T0 = np.zeros(nc)
            T0[j] = 1.0
            for i in range(1,nround+1):
                T0 = np.dot(T,T0)
                TPs[j,i] = T0[j]
        outf = os.path.join(self.outp, "ProbMSM")
        self.MSMTP = TPs
        self._flag = True

    def plotv1(self,md,msm,var,title,name):
        fig,ax = plt.subplots()
        if self.unit=='ps':
            x = [self.lagtime*i*self.dt*1.0 for i in range(self.nround+1)]
        elif self.unit=='ns':
            x = [self.lagtime*i*self.dt/1000.0 for i in range(self.nround+1)]
        ax.errorbar(x, md, yerr=var,fmt="o",label='MD',markerfacecolor='none', elinewidth=3.0, c="red",capsize=5)
        ax.scatter(x,md, color=None, edgecolors="red",s=30)
        ax.plot(x, msm, '-', label='MSM',lw=4.0, c="blue")
        ax.set_xlabel('lag time (%s)'%self.unit)
        ax.set_ylabel('residence probability')
        ax.set_title(title)
        #plt.text(x[-1]/2,0.8,title)
        ax.set_ylim((0,1))
        ax.set_xlim(0,x[-1]+self.lagtime*self.dt/4000)
        fig.tight_layout()
        fig.savefig(os.path.join(self.outp,name),dpi=100)
        plt.close()

    def plotv2(self,md,msm,var,title,ax):
        if self.unit=='ps':
            lag = self.lagtime*self.dt*1.0
        elif self.unit=='ns':
            lag = self.lagtime*self.dt/1000.0
        else:
            print("ERROR: not supported units: %s"%self.unit)
        x = [lag*i for i in range(self.nround+1)]
        ax.errorbar(x, md, yerr=var, fmt='none', label='MD', elinewidth=2,markersize=10,capsize=10,mfc=None,ecolor="red",color="red")
        #ax.scatter(x, md)
        ax.scatter(x,md,c="none",edgecolors="red",s=50,linewidths =1.5)
        ax.plot(x, msm, '-', label='MSM',lw=3.0, c="blue")
        ymax = max([max(md),max(msm)])
        ax.text(x[-1]/2,1.15,title, verticalalignment='top',horizontalalignment='center')

        ax.set_xlim(0-lag/4*3,x[-1]+lag/2)
        if len(x)<=4:
            ax.set_xticks(x)
        else:
            n = len(x)
            ax.set_xticks(x[:n:3])
            #x = [0,x[n-1]/3,x[n-1]/3*2,x[n-1]]
            #ax.set_xticks(x)
        ax.set_yticks(np.arange(0,ymax+0.5,0.5))

    def Cktest_plot(self,P,top=10):
        P = [(i,p) for i,p in enumerate(P)]
        self.P = self.pi

        P = sorted(P,key=lambda a:a[-1],reverse=True)
        if not self._flag:
            raise Exception("Please do the CKtest first.")
        #for i in range(top):
        #    md = self.MDTP[P[i][0],:]
        #    msm = self.MSMTP[P[i][0],:]
        #    var = self.VAR[P[i][0],:]
        #    title = r'$state=%d,sigma=%.4f$'%tuple(P[i])
        #    name = "%d_%d.pdf"%(i,P[i][0])
        #    self.plotv1(md,msm,var,title,name)
        '''
        if top<10:
            for i in range(len(P)):
                md = self.MDTP[P[i][0],:]
                msm = self.MSMTP[P[i][0],:]
                var = self.VAR[P[i][0],:]
                title = r'P=%.3f'%P[i][1]
                name = "%d_%d.pdf"%(i,P[i][0])
                self.plotv1(md,msm,var,title,name)
            return 0
        '''
        fig,axs = plt.subplots(3,3,sharex=True,sharey=True)#,constrained_layout=True)
        fig.set_size_inches(6,5)
        axs = axs.flatten()
        for i,ax in enumerate(axs):
            if i<top:
                md = self.MDTP[P[i][0],:self.nround+1]
                msm = self.MSMTP[P[i][0],:self.nround+1]
                var = self.VAR[P[i][0],:self.nround+1]
                #title = u'$S%d,p%.3f$'%tuple(P[i])
                title = u'P=%.3f'%P[i][1]
                self.plotv2(md,msm,var,title,ax)
            if i == 3:
                ax.set_ylabel('residence probability')
            if i == 7:
                ax.set_xlabel('lag time (%s)'%self.unit)
        ax.set_ylim(0, 1.2)
        #fig.set_constrained_layout_pads(h_pad=2./72.,w_pad=2./72.,wspace=0.,hspace=0.)
        fig.subplots_adjust(left=0.15,bottom=0.15,top=0.9,right=0.95,hspace=0.2,wspace=0.25)
        #fig.subplots_adjust(left=0.,bottom=0.,top=0.1,right=0.1)
        #fig.tight_layout()
        plt.show()
        fig.savefig(os.path.join(self.outp,'combine.png'), dpi=300)

    def save(self):
        if not self._flag:
            raise Exception("Please do the CKtest first.")

        np.savetxt(self.MDTPf,self.MDTP)

        np.savetxt(self.MSMTPf,self.MSMTP)

        np.savetxt(self.VARf,self.VAR)
        Pf = os.path.join(self.outp,'Populations')

        np.savetxt(Pf, self.P)

def parse_args():

    parser = argparse.ArgumentParser(description="chapman-kolmogorov test for msm")
    parser.add_argument("-m", dest="msmf", help="The msm pickl file, from msmbuilder 3.*. or from pyemma", required=True)
    parser.add_argument("-s", dest="seq", help="The cluster label pickle file. required for msmbuilder msmf input.", default=None)
    parser.add_argument("-dt", help="The time of every frame. default=200 ps", default=200, type=int)
    parser.add_argument("-du", help="The units of the dt (ps, ns, ms...) , default: ps", default="ps")
    parser.add_argument("-nr", help="The round to cal CK-test. default=3", default=3, type=int)
    parser.add_argument("-o", dest="outp", help="The result path to save output. default: cktest", default="cktest")
    parser.add_argument("-T", help="The number of thread to run. default: 4", type=int, default=4)

    args = parser.parse_args()

    return args.msmf,  args.dt, args.du, args.nr, args.outp, args.T, args.seq


def convert_k_v(dic):
    new_dic = {}
    for k,v in dic.items():
        new_dic[v] = k

    return new_dic

def main():
    msmf, dt, du, nr, outp, T, seqf = parse_args()

    # load pickl data
    msm = pd.read_pickle(msmf)
    if hasattr(msm,'n_states_'):
        if not seqf:
            print("ERROR: msmbuilder msm pickle requires the cluster assign file.")
            sys.exit(0)
        seq_original = pd.read_pickle(seqf)
        msm.mapping_ = convert_k_v(msm.mapping_)
        seq = msm.transform(seq_original, mode="fill")
        #print(set(np.concatenate(seq)))
        nstates,Pi, lagtime = msm.n_states_,msm.populations_, msm.lag_time
        matrix = msm.transmat_
    else:
        seq = msm.discrete_trajectories_full
        nstates,Pi,lagtime = msm.nstates, msm.pi, msm.lagtime
        matrix = msm.transition_matrix

        maxstates = np.concatenate(seq).max()+1
        if len(msm.active_set)!=maxstates:
            mapping = {}
            k = 0
            FLAG = False
            for i in range(maxstates):
                if i not in msm.active_set:
                    mapping[i] = np.nan 
                    FLAG = True
                else:
                    mapping[i] = k
                    k += 1
            if FLAG:
                seq_new = []
                for s in seq:
                    s = s.astype("float")
                    for k,v in mapping.items():
                        mask = s==k
                        #print mask
                        s[mask] = v 
                        seq_new.append(s)
                seq = np.array(seq_new)
    ck = ChapmanKolmogorovTest(msm, dt, du, outp, T, nstates, lagtime, Pi, nr)
    ck.CK_MSM(matrix)
    ck.CK_MD(seq)
    # plot
    if nstates >50:
        n = 50
    else:
        n = nstates
    ck.Cktest_plot(Pi, top=n)
    # save
    ck.save()

if __name__ == "__main__":
    main()
