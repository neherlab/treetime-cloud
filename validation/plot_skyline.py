from external_binaries import FFPOPSIM_SKYLINE_BIN
from utility_functions_simulated_data import _ffpopsim_tree_aln_postprocess, generations_from_ffpopsim_tree
import sys,os, glob

import numpy as np
from treetime import TreeTime
from matplotlib import pyplot as plt
import seaborn as sns

from plot_defaults import *
cols = sns.color_palette(n_colors=3)

def read_estimate_skyline(infile):
    with open(infile, 'r') as inf:
        ss = inf.readlines()

    period = None
    amp = None
    x,y,trueY = [],[],[]

    reading_arrays = False
    for s in ss:
        if s.startswith('#period'):
            period = float(s.split('=')[-1].strip())
        elif s.startswith('#amplitude'):
            amp = float(s.split('=')[-1].strip())
        elif s.startswith('#x, y, trueY'):
            reading_arrays = True
        elif reading_arrays:
            try:
                i,j,k = map(float, s.split(','))
                x.append(i)
                y.append(j)
                trueY.append(k)
            except:
                reading_arrays = False
        else:
            continue
    return period, amp, x, y, trueY

def plot_skyline(res, periods, bottle_necks, figname=None):
    fig = plt.figure(figsize=onecolumn_figsize)

    # common Y label
    fig.text(0.05,0.5, 'Coalescent population size estimate',
        rotation='vertical',
        fontsize=label_fs,
        verticalalignment='center',
        horizontalalignment='center')

    axs = [fig.add_subplot(211), fig.add_subplot(212)]
    for ri, (p, amp, x, s, t) in enumerate(res):

        if not p in periods:
            continue

        axidx = periods.index(p)
        colidx = bottle_necks.index(amp)
        ax = axs[axidx]
        print(p,amp, np.corrcoef(s,t)[0,1])
        ax.plot(x, s, c=cols[colidx], ls='-', label ="bottleneck: %1.1f"%(((1-amp))))
        ax.plot(x, t, c=cols[colidx], ls='--')
        ax.set_xlim(4300,4800)


        for label in ax.get_yticklabels():
            label.set_fontsize(tick_fs)

        if axidx == 0:
            for label in ax.get_xticklabels():
                label.set_visible(False)
        else:
            for label in ax.get_xticklabels():
                label.set_fontsize(tick_fs)
    axs[0].legend(loc=2, fontsize=legend_fs)

    #axes.set_ylabel('Population size estimate', fontsize=label_fs)
    axs[1].set_xlabel('Time in simulated generations', fontsize = label_fs)

    if figname is not None:
        for fmt in formats:
            plt.savefig('{}.{}'.format(figname, fmt))


if __name__ == '__main__':

    SAVE_FIG = True

    resdir = './skyline/'

    periods = [0.5, 2.0]
    bottle_necks = [0.5, 0.8, 0.9]

    fnames = glob.glob(os.path.join(resdir,'FF*Mu0.001*res.txt'))
    res = map(read_estimate_skyline, fnames)

    plot_skyline(res, periods, bottle_necks,
        figname='./figs/skyline' if SAVE_FIG else None)


