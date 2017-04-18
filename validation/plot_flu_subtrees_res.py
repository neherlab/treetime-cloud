import pandas
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import os, sys
import shutil
from Bio import Phylo

import utility_functions_flu as flu_utils
import utility_functions_beast as beast_utils

from plot_defaults import *

##  Read datasets as-is
def read_lsd_dataset(fname):
    """
    TODO
    """
    lsd_cols = ['File', 'N', 'Tmrca_sim', 'mu_sim', 'Runtime', 'objective']
    lsd_df = pandas.read_csv(fname, names=lsd_cols)
    return lsd_df

def read_treetime_dataset(fname):
    """
    TODO
    """
    cols = ['File', 'N', "Tmrca_sim", "mu_sim", "R2_leaves", "R2_internal", "Runtime"]
    df = pandas.read_csv(fname, names=cols)
    return df

##  Make statistics on the pure datasets
def beast_logs_stat(logsdir, treedir):
    """
    logsdir: directory with beast log files
    treedir: directory with the initial trees used for the beast run
    """

    Ns = []
    LH_mean = []
    LH_err = []
    Tmrca_mean = []
    Tmrca_err = []
    Mu_mean = []
    Mu_err = []

    logsfiles = [k for k in os.listdir(logsdir) if k.endswith('.log.xml')]

    for log in logsfiles:

        treename = os.path.split(log)[-1][:-8] + '.nwk'
        treepath = os.path.join(treedir, treename)
        logpath = os.path.join(logsdir, log)
        dates = flu_utils.dates_from_flu_tree(treepath)

        df = beast_utils.read_beast_log(logpath, np.max(dates.values()))
        if df is None or df.shape[0] < 300:
            continue
        new_Ns = Phylo.read(treepath, 'newick').count_terminals()
        if new_Ns > 1e3:
            continue
        Ns.append(new_Ns)

        LH_mean.append(df['likelihood'].mean())
        LH_err.append(df['likelihood'].std())

        Mu_mean.append(df['clock.rate'].mean())
        Mu_err.append(df['clock.rate'].std()) # * Mu_mean[-1] / 2.)

        Tmrca_mean.append(df['treeModel.rootHeight'].mean())
        Tmrca_err.append(df['treeModel.rootHeight'].std())

    res = pandas.DataFrame({
        "Ns" : Ns,
        "Tmrca_mean" : Tmrca_mean,
        "Tmrca_err" : Tmrca_err,
        "Mu_mean" : Mu_mean,
        "Mu_err" : Mu_err,
        "LH_mean" : LH_mean,
        "LH_err" : LH_err
        })

    res = res.sort('Ns')
    return res

def treetime_dataset_stat(df):

    Tmrca_mean = []
    Tmrca_err = []
    Mu_mean = []
    Mu_err = []
    Runtime_mean = []
    Runtime_err = []

    Ns = df["N"].unique()

    Nsidx = np.ones(Ns.shape, dtype=bool)
    for idx, N in enumerate(Ns):
        Nidx = df["N"] == N
        if Nidx.sum() == 0:
            Nsidx[idx] = False
            continue
        Tmrca_mean.append(df[Nidx]["Tmrca_sim"].mean())
        Tmrca_err.append(df[Nidx]["Tmrca_sim"].std())
        Mu_mean.append(df[Nidx]["mu_sim"].mean())
        Mu_err.append(df[Nidx]["mu_sim"].std())
        Runtime_mean.append(df[Nidx]["Runtime"].mean())
        Runtime_err .append(df[Nidx]["Runtime"].std())

    df = pandas.DataFrame({
        "Ns" : Ns[Nsidx],
        "Tmrca_mean" : Tmrca_mean,
        "Tmrca_err" : Tmrca_err,
        "Mu_mean" : Mu_mean,
        "Mu_err" : Mu_err,
        "Runtime_mean" : Runtime_mean,
        "Runtime_err" : Runtime_err
        })
    return df

def lsd_dataset_stat(df):

    Tmrca_mean = []
    Tmrca_err = []
    Mu_mean = []
    Mu_err = []
    Runtime_mean = []
    Runtime_err = []

    Ns = df["N"].unique()
    Nsidx = np.ones(Ns.shape, dtype=bool)

    for idx, N in enumerate(Ns):
        Nidx = df["N"] == N
        if Nidx.sum() == 0:
            Nsidx[idx] = False
            continue

        Tmrca_mean.append(df[Nidx]["Tmrca_sim"].mean())
        Tmrca_err.append(df[Nidx]["Tmrca_sim"].std())
        Mu_mean.append(df[Nidx]["mu_sim"].mean())
        Mu_err.append(df[Nidx]["mu_sim"].std())
        Runtime_mean.append(df[Nidx]["Runtime"].mean())
        Runtime_err .append(df[Nidx]["Runtime"].std())

    res = pandas.DataFrame({
        "Ns" : Ns[Nsidx],
        "Tmrca_mean" : Tmrca_mean,
        "Tmrca_err" : Tmrca_err,
        "Mu_mean" : Mu_mean,
        "Mu_err" : Mu_err,
        "Runtime_mean" : Runtime_mean,
        "Runtime_err" : Runtime_err
        })

    return res

## Plot statistics
def plot_res(what, tt=None, lsd=None, beast=None, save=True, suffix=None, **kwargs):

    if what == 'Tmrca':
        mean = 'Tmrca_mean'
        err = 'Tmrca_err'
        title = "Estimated Tmrca as function of sample size\nLSD params: -{}".format(suffix)

        ylabel = "$T_{mrca}, [year]$"

    elif what == "Mu":
        mean = 'Mu_mean'
        err =  'Mu_err'
        title =  "Estimated Mutation rate as function of sample size\nLSD params: -{}".format(suffix)

        ylabel = "Mutation rate, $[\mathrm{year}^{-1}]$"

    fig = plt.figure(figsize=onecolumn_figsize)
    axes = fig.add_subplot(111)
    axes.ticklabel_format(useOffset=False)

    if tt is not None:
        axes.errorbar(tt['Ns'], tt[mean], tt[err]/2, markersize=markersize, marker='o', c='b', label='TreeTime')

    if lsd is not None:
        axes.errorbar(lsd['Ns'], lsd[mean], lsd[err]/2, markersize=markersize, marker='o', c='g', label='LSD')

    if beast is not None:
        axes.errorbar(beast['Ns'], beast[mean], beast[err]/2, markersize=markersize, marker='o', c='r', label='BEAST')

    axes.grid('on')
    axes.legend(loc=0)
    axes.set_xscale('log')

    axes.set_ylabel(ylabel, fontsize=label_fs)
    axes.set_xlabel("Tree size, #sequences",fontsize=label_fs)
    #axes.set_title(title)

    for label in axes.get_xticklabels(): 
            label.set_fontsize(tick_fs) 
    for label in axes.get_yticklabels(): 
            label.set_fontsize(tick_fs)

    if save:
        fig.savefig("./figs/{}_LSD_{}.svg".format(what, suffix))
        fig.savefig("./figs/{}_LSD_{}.jpg".format(what, suffix))

if __name__ == "__main__":

    #  directory to search for the result tables:
    res_dir = './flu_H3N2/subtree_samples/'
    treetime_res = os.path.join(res_dir, 'treetime_res.csv')
    lsd_res = os.path.join(res_dir, 'lsd_res.csv')
    beast_log_dir = os.path.join(res_dir, 'beast_out_cluster_2')
    beast_tree_dir = os.path.join(res_dir, 'subtrees')


    PLOT_TREETIME = True
    PLOT_LSD = False
    PLOT_BEAST = True


    if PLOT_TREETIME:
        tt_df = treetime_dataset_stat(read_treetime_dataset(treetime_res))
    else:
        tt_df = None
    if PLOT_LSD:
        lsd_df = lsd_dataset_stat(read_lsd_dataset(lsd_res))
    else:
        lsd_df = None
    if PLOT_BEAST:
        beast = beast_logs_stat(beast_log_dir, beast_tree_dir)
    else:
        beast=None

    plot_res('Tmrca', tt=tt_df, lsd=lsd_df, beast=beast, save=True)
    plot_res('Mu', tt=tt_df, lsd=lsd_df, beast=beast, save=True)

