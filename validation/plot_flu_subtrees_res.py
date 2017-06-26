import pandas
import numpy as np
#from scipy.stats import linregress
import matplotlib.pyplot as plt
import os
#import shutil
#from Bio import Phylo

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

    logsfiles = [k for k in os.listdir(logsdir) if k.endswith('.log.txt')]

    for log in logsfiles:

        treename = os.path.split(log)[-1][:-8] + '.nwk'
        treepath = os.path.join(treedir, treename)
        logpath = os.path.join(logsdir, log)
        dates = flu_utils.dates_from_flu_tree(treepath)

        df = beast_utils.read_beast_log(logpath, np.max(dates.values()))
        if df is None or df.shape[0] < 300:
            continue
        new_Ns = int(treename.split('_')[-2])
        #new_Ns = Phylo.read(treepath, 'newick').count_terminals()
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

def make_beast_pivot(beast_res):
    out_Ns = np.unique(beast_res['Ns'])
    out = pandas.DataFrame({
        "Ns" : out_Ns,
        "Tmrca_mean" : [np.mean(beast_res[beast_res['Ns'] == x]['Tmrca_mean'])
                            for x in out_Ns],
        "Tmrca_err" : [np.std(beast_res[beast_res['Ns'] == x]['Tmrca_mean'])
                            for x in out_Ns],
        "Mu_mean" : [np.mean(beast_res[beast_res['Ns'] == x]['Mu_mean'])
                            for x in out_Ns],
        "Mu_err" : [np.std(beast_res[beast_res['Ns'] == x]['Mu_mean'])
                            for x in out_Ns],
        "LH_mean" : [np.mean(beast_res[beast_res['Ns'] == x]['LH_mean'])
                            for x in out_Ns],
        "LH_err" : [np.std(beast_res[beast_res['Ns'] == x]['LH_mean'])
                            for x in out_Ns],
        })
    return out

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
def plot_res(what, tt=None, lsd=None, beast=None, save=True, suffix=None, scatter_points=True, **kwargs):

    if what == 'Tmrca':
        mean = 'Tmrca_mean'
        err = 'Tmrca_err'
        #title = "Estimated Tmrca as function of sample size\nLSD params: -{}".format(suffix)

        ylabel = "T$\mathrm{_{mrca}}, [\mathrm{Year}]$"

    elif what == "Mu":
        mean = 'Mu_mean'
        err =  'Mu_err'
        #title =  "Estimated substitution rate as function of sample size\nLSD params: -{}".format(suffix)

        ylabel = "Substitution rate, [$\mathrm{Year}^{-1}$]"

    fig = plt.figure(figsize=onecolumn_figsize)
    axes = fig.add_subplot(111)
    axes.ticklabel_format(useOffset=False)
    axes.set_xscale('log')

    if tt is not None:
        if scatter_points:
            x, y = shift_point_by_markersize (axes, tt['Ns'], tt[mean], markersize/2.0)
        else:
            x, y = tt['Ns'], tt[mean]
        axes.errorbar(x, y, tt[err]/2, markersize=markersize, marker='o', c=tt_color, label='TreeTime')

    if lsd is not None:
        if scatter_points:
            x, y = shift_point_by_markersize (axes, lsd['Ns'], lsd[mean], -1.*markersize/2.0)
        else:
            x, y = lsd['Ns'], lsd[mean]
        axes.errorbar(x, y, lsd[err]/2, markersize=markersize, marker='o', c=lsd_color, label='LSD')

    if beast is not None:
        #  beast points stay in the center
        x, y = beast['Ns'], beast[mean]
        # if scatter_points:
        #     shift_point_by_markersize (axes, beast['Ns'], beast[mean], -1.*markersize/2.0)
        # else:
        #     x, y = beast['Ns'], beast[mean]
        axes.errorbar(x, y, beast[err]/2, markersize=markersize, marker='o', c=beast_color, label='BEAST')

    axes.grid('on')
    axes.legend(loc=0,fontsize=legend_fs)
    axes.set_ylabel(ylabel, fontsize=label_fs)
    axes.set_xlabel("Tree size, [#$\mathrm{Sequences}$]",fontsize=label_fs)
    #axes.set_title(title)

    for label in axes.get_xticklabels():
            label.set_fontsize(tick_fs)
    for label in axes.get_yticklabels():
            label.set_fontsize(tick_fs)

    if save:
        fig.savefig("./figs/fluH3N2_subtrees_{}.svg".format(what))
        fig.savefig("./figs/fluH3N2_subtrees_{}.png".format(what))
        fig.savefig("./figs/fluH3N2_subtrees_{}.pdf".format(what))

if __name__ == "__main__":

    #  directory to search for the result tables:
    res_dir = './flu_H3N2/subtree_samples/'
    treetime_res = os.path.join(res_dir, '2017-04-20_treetime_res.csv')
    lsd_res = os.path.join(res_dir, '2017-04-20_lsd_res.csv')
    beast_log_dir = os.path.join(res_dir, '2017-04-20/beast_out')
    beast_tree_dir = os.path.join(res_dir, '2017-04-20/subtrees')

    PLOT_TREETIME = True
    PLOT_LSD = True
    PLOT_BEAST = True

    if PLOT_TREETIME:
        tt_df = treetime_dataset_stat(read_treetime_dataset(treetime_res))
        tt_df = tt_df.sort(columns='Ns')
    else:
        tt_df = None
    if PLOT_LSD:
        lsd_df = lsd_dataset_stat(read_lsd_dataset(lsd_res))
        lsd_df = lsd_df.sort(columns='Ns')
    else:
        lsd_df = None
    if PLOT_BEAST:
        beast = beast_logs_stat(beast_log_dir, beast_tree_dir)
        beast = make_beast_pivot(beast)
    else:
        beast=None

    plot_res('Tmrca', tt=tt_df, lsd=lsd_df, beast=beast, save=True)
    plot_res('Mu', tt=tt_df, lsd=lsd_df, beast=beast, save=True)

