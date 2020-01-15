import numpy as np
from Bio import AlignIO, Phylo
from Bio.Align import  MultipleSeqAlignment
import random
import subprocess
import datetime
import os, copy
import matplotlib.pyplot as plt
from scipy.stats import linregress
from collections import Counter
import treetime
import pandas
from plot_defaults import *
from scipy.interpolate import interp1d
from plot_defaults import shift_point_by_markersize
import utility_functions_beast as beast_utils
import utility_functions_flu as flu_utils

def read_treetime_csv(inf):
    cols = ['File', 'Frac', 'Tmrca', 'Mu', 'Mu_R2', 'Internal_R2', 'Runtime']
    df = pandas.read_csv(inf, names=cols, header=0)
    return df

def read_beast_csv(fname):
    cols = ['Filename','Frac','LH','LH_std','Tmrca','Tmrca_std','Mu','Mu_std']
    df = pandas.read_csv(fname, names=cols, header=0)
    return df

def make_results_pivot(df):

    fs = np.unique(df['Frac'])
    fs = fs[fs >= 0.1]
    Tmrca = []
    Tmrca_err = []
    Mu = []
    Mu_err = []
    R2 = []
    R2_err = []

    for frac in fs:
        idxs = df['Frac'] == frac
        if idxs.sum() < 10:
            continue

        Tmrca.append((df['Tmrca'][idxs]).mean())
        Tmrca_err.append((df['Tmrca'][idxs]).std())
        Mu.append((df['Mu'][idxs]).mean())
        Mu_err.append((df['Mu'][idxs]).std())
        R2.append((df['Mu_R2'][idxs]).mean())
        R2_err.append((df['Mu_R2'][idxs]).std())

    res = pandas.DataFrame(
        {'Frac': fs,
         'Tmrca': Tmrca,
         'Tmrca_err': Tmrca_err,
         'Mu': Mu,
         'Mu_err': Mu_err,
         'R2': R2,
         'R2_err': R2_err})
    return res

def make_beast_pivot(beast_res):

    out_Ns = np.unique(beast_res['Frac'])
    out = pandas.DataFrame({
        "Frac" : out_Ns,
        "Tmrca" : [np.mean(beast_res[beast_res['Frac'] == x]['Tmrca'])
                            for x in out_Ns],
        "Tmrca_err" : [np.std(beast_res[beast_res['Frac'] == x]['Tmrca'])
                            for x in out_Ns],
        "Mu" : [np.mean(beast_res[beast_res['Frac'] == x]['Mu'])
                            for x in out_Ns],
        "Mu_err" : [np.std(beast_res[beast_res['Frac'] == x]['Mu'])
                            for x in out_Ns],
        "LH" : [np.mean(beast_res[beast_res['Frac'] == x]['LH'])
                            for x in out_Ns],
        "LH_err" : [np.std(beast_res[beast_res['Frac'] == x]['LH'])
                            for x in out_Ns],
        })
    return out

def plot_results(what, treetime=None, beast=None, figname=None):

    if what == 'Tmrca':
        ylabel = r'$\mathrm{T_{mrca}}, [\mathrm{Year}]$'
        title = 'Dependence of $\mathrm{T}_{mrca}$ prediction on the fraction of the known dates'

    elif what == 'Mu':
        ylabel = 'Mutation rate, $[\mathrm{Year}^{-1}]$'
        title = 'Dependence of Mutation rate assessment on the fraction of the known dates'

    elif what == 'R2':
        ylabel = 'Molecular clock regression'
        title = 'Quality of the molecular clock regression vs fraction of known dates'


    fig = plt.figure(figsize=onecolumn_figsize)
    axes = fig.add_subplot(111)

    if treetime is not None:
        x, y = shift_point_by_markersize(axes, treetime.Frac, treetime[what], -markersize/2)
        axes.errorbar(x, y, treetime[what + '_err'], fmt='o-', label="TreeTime", markersize=markersize, c=tt_color)


    if beast is not None:
        x, y = shift_point_by_markersize(axes, beast.Frac, beast[what], markersize/2)
        axes.errorbar(x, y, beast[what + '_err'], fmt='o-', label="BEAST", markersize=markersize, c=beast_color)

    axes.ticklabel_format(useOffset=False)
    axes.set_xlabel('Fraction of known dates', fontsize = label_fs)
    axes.set_ylabel(ylabel, fontsize=label_fs)
    #axes.set_title(title)
    axes.grid('on')
    axes.legend(loc=0, fontsize=legend_fs)
    for label in axes.get_xticklabels():
            label.set_fontsize(tick_fs)
    for label in axes.get_yticklabels():
            label.set_fontsize(tick_fs)

    if figname is not None:
        for fmt in formats:
            fig.savefig("{}.{}".format(figname, fmt))

if __name__ == '__main__':

    SAVE_FIG = True

    ##
    ## CSV tables with the results:
    ##
    beast_csv = "./flu_H3N2/missing_dates/beast_res.csv"
    treetime_csv = "./flu_H3N2/missing_dates/treetime_res.csv"

    ##
    ## Read results, create pivot tables:
    ##
    beast_pivot = make_beast_pivot(read_beast_csv(beast_csv))
    treetime_pivot = make_results_pivot(read_treetime_csv(treetime_csv))

    ##
    ## Plot the results:
    ##
    plot_results('Tmrca', treetime=treetime_pivot, beast=beast_pivot,
        figname="./figs/fluH3N2_missingDates_100seqs" if SAVE_FIG else None)

    plot_results('Mu', treetime=treetime_pivot, beast=beast_pivot,
        figname="./figs/fluH3N2_missingDates_100seqs" if SAVE_FIG else None)
