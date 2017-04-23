import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import os, sys
import pandas
from Bio import Phylo

import utility_functions_beast as beast_utils
import utility_functions_simulated_data as sim_utils

from plot_defaults import *

def read_treetime_results_dataset(fname):
    """
    Read results of the TreeTime simulations

    Args:
     - fname: path to the input file

    Returns:
     - df: Table of results as pandas data-frame
    """

    columns = ['File', 'Sim_Tmrca', 'Tmrca', 'mu', 'R', 'R2_int']
    df = pandas.read_csv(fname, names=columns)

    #filter obviously failed simulations
    df = df[[len(str(k)) > 10 for k in df.File]]
    df = df[df.R > 0.1]

    # some very basic preprocessing
    df['dTmrca'] = (df['Sim_Tmrca'] - df['Tmrca'])
    df['Sim_mu'] = map(lambda x: float(x.split("/")[-1].split('_')[6][2:]), df.File)
    df['Ns'] = map(lambda x: int(x.split("/")[-1].split('_')[3][2:]), df.File)
    df['Ts'] = map(lambda x: int(x.split("/")[-1].split('_')[4][2:]), df.File)
    df['N'] = map(lambda x: int(x.split("/")[-1].split('_')[2][1:]), df.File)
    df['T'] = df['Ns']*df['Ts']
    df['Nmu'] = (df['N']*df['Sim_mu'])

    return df

def read_lsd_results_dataset(fname):
    """
    Read results of the LSd simulations

    Args:
     - fname: path to the input file

    Returns:
     - df: Table of results as pandas data-frame
    """

    columns = ['File', 'Sim_Tmrca', 'Tmrca', 'mu', 'obj']
    df = pandas.read_csv(fname, names=columns)

    # Filter out obviously wrong data
    df = df[[len(k) > 10 for k in df.File]]

    #Some basic preprocessing
    df['dTmrca'] = (df['Sim_Tmrca'] - df['Tmrca'])
    df['Sim_mu'] = map(lambda x: float(x.split("/")[-1].split('_')[6][2:]), df.File)
    df['Ns'] = map(lambda x: int(x.split("/")[-1].split('_')[3][2:]), df.File)
    df['Ts'] = map(lambda x: int(x.split("/")[-1].split('_')[4][2:]), df.File)
    df['N'] = map(lambda x: int(x.split("/")[-1].split('_')[2][1:]), df.File)
    df['T'] = df['Ns']*df['Ts']
    df['Nmu'] = (df['N']*df['Sim_mu'])

    return df

def create_lsd_tt_pivot(df, T_over_N=None, mean_or_median='median'):

    if T_over_N is not None:
        DF = df[ df["T"] / df["N"] == T_over_N ]
    else:
        DF = df

    N_MUS = np.unique(DF.Nmu)
    N_MUS_idxs = np.ones(N_MUS.shape, dtype=bool)

    mu_mean = []
    mu_err = []
    tmrca_mean = []
    tmrca_err = []

    for idx, N_MU  in enumerate(N_MUS):
        idxs = DF.Nmu == N_MU
        if idxs.sum() == 0:
            N_MUS_idxs[idx] = False
            continue

        dMu = (DF.Sim_mu[idxs] - DF.mu[idxs])/DF.Sim_mu[idxs]
        dMu.sort_values(inplace=True)
        dMu[int(dMu.shape[0]*0.05) : int(dMu.shape[0]*0.95)]

        dTmrca = DF.dTmrca[idxs]/DF.N[idxs]
        dTmrca.sort_values(inplace=True)
        dTmrca[int(dTmrca.shape[0]*0.05) : int(dTmrca.shape[0]*0.95)]

        if mean_or_median == "mean":
            mu_mean.append(np.mean(dMu))
            mu_err.append(np.std(dMu))

            tmrca_mean.append(np.mean(dTmrca))
            tmrca_err.append(np.std(dTmrca))
        else:
            q75, q25 = np.percentile(dMu, [75 ,25])
            mu_err.append((q75 - q25)) #np.std(DF.dTmrca[idxs])
            mu_mean.append(np.median(dMu))
            q75, q25 = np.percentile(dTmrca, [75 ,25])
            tmrca_err.append((q75 - q25)) #np.std(DF.dTmrca[idxs])
            tmrca_mean.append(np.median(dTmrca))


    res = pandas.DataFrame({
        "Nmu" : N_MUS[N_MUS_idxs],
        "dMu_mean" : mu_mean,
        "dMu_err" : mu_err,
        "dTmrca_mean" : tmrca_mean,
        "dTmrca_err" : tmrca_err,
        })
    res = res.sort_values(by='Nmu')
    return res

def plot_raw_data(df):
    """
    Plot the raw simulations data, clustered by the $N\mu$ value
    """

    # number of mutation classes
    N_MUS = np.unique(df.Nmu)

    #  set colors
    NUM_COLORS = len(N_MUS)
    cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)

    # plot results
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(NUM_COLORS)])

    for MU in N_MUS:
        #ax.plot(df["T"]/df['N'], df.dTmrca / (2016.5 - df.Tmrca), 'o', label="MU = " + str(MU), alpha=0.7)
        ax.plot(df["T"]/df['N']*(1+np.random.normal(size=len(df))*0.1), df.dTmrca / df["N"], 'o', label="MU = " + str(MU), alpha=0.6)

    ax.legend(loc=0)
    plt.xlabel("Total evolution time,  $T/N$")
    plt.ylabel("Relative error in Tmrca estimation, $\Delta T_{mrca} / N$")

def plot_mutation_rate_distributions(df, TN_min=3, TN_max=10, plot_title=""):
    """
    TODO
    """

    # number of mutation classes
    N_MUS = np.unique(df.Nmu)


    #  set colors
    NUM_COLORS = len(N_MUS)
    cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)

    DF = df[ (df['T']/df['N'] > TN_min) & (df['T']/df['N'] < TN_max) ]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(NUM_COLORS)])




    for idx,MU  in enumerate(N_MUS):
        idxs = DF.Nmu == MU
        print (MU, np.sum(idxs))
        if idxs.sum() == 0:
            continue

        # plot histogram
        dMu = (DF.Sim_mu[idxs] - DF.mu[idxs])/DF.Sim_mu[idxs]
        plt.hist(dMu, label="$N\cdot\mu = $" + str(MU), alpha=.5, bins=20)

    ax.set_title(plot_title)
    plt.legend(loc=0)

def T_over_N_from_filename(filename):
    return 1. * int(filename.split('_')[4][2:]) * int(filename.split('_')[3][2:]) / int(filename.split('_')[2][1:])

def read_beast_log(filename):
    return beast_utils.read_beast_log(filename, sim_utils.NEAREST_DATE)

def read_all_beast_logs(logsdir, treesdir, T_over_N=None):

    def Tmrca_from_tree(logfile, treesdir):
        treename = os.path.join(treesdir, logfile[:-8] + ".nwk")
        Tmrca, dates = sim_utils.dates_from_ffpopsim_tree(treename)
        return Tmrca

    if T_over_N is not None:
        beast_logs = [k for k in os.listdir(logsdir) if k.endswith('.log.txt')
                                     and T_over_N_from_filename(k) == T_over_N]
    else:
        beast_logs = [k for k in os.listdir(logsdir) if k.endswith('.log.txt')]

    File = []
    Sim_Tmrca = []
    dTmrca = []
    Sim_Mu = []
    Ns = []
    Ts = []
    N = []
    T = []
    Nmu = []

    LH = []
    LH_std = []
    Tmrca = []
    Tmrca_std = []
    Mu = []
    Mu_std = []
    dMu = []

    for beast_log in beast_logs:
        df = read_beast_log(os.path.join(logsdir, beast_log))
        if df is None or df.shape[0] < 200 :
            print ("Beast log {} is BAD".format(beast_log))
            continue

        File.append(beast_log)

        Sim_Tmrca.append(Tmrca_from_tree(beast_log, treesdir))
        Sim_Mu.append(float(beast_log.split("/")[-1].split('_')[6][2:]))
        Ns.append(int(beast_log.split("/")[-1].split('_')[3][2:]))
        Ts.append(int(beast_log.split("/")[-1].split('_')[4][2:]))
        N.append(int(beast_log.split("/")[-1].split('_')[2][1:]))
        T.append(Ns[-1] * Ts[-1])
        Nmu.append(N[-1]*Sim_Mu[-1])

        LH.append(df['likelihood'][-50:].mean())
        LH_std.append(df['likelihood'][-50:].std())
        Tmrca.append(df['treeModel.rootHeight'][-50:].mean())
        Tmrca_std.append(df['treeModel.rootHeight'][-50:].std())
        Mu.append(df['clock.rate'][-50:].mean())
        Mu_std.append(df['clock.rate'][-50:].std())

        dTmrca.append(Sim_Tmrca[-1] - Tmrca[-1])
        dMu.append(Sim_Mu[-1] - Mu[-1])

    res = pandas.DataFrame({
        "File" : File,
        "Likelihood" : LH,
        "Likelihood_std" : LH_std,
        "Tmrca" : Tmrca,
        "Tmrca_std" : Tmrca_std,
        "Mu": Mu,
        "Mu_std" : Mu_std,
        "Sim_Tmrca" : Sim_Tmrca,
        "Sim_Mu" : Sim_Mu,
        "Ns" : Ns,
        "Ts" : Ts,
        "N" : N,
        "T" : T,
        "Nmu" : Nmu,
        "dTmrca" : dTmrca,
        "dMu" : dMu

        })

    return res

def create_beast_log_pivot(df, T_over_N=None):


    if T_over_N is not None:
        DF = df[ df["T"] / df["N"] == T_over_N ]
    else:
        DF = df

    mu_mean = []
    mu_err = []
    tmrca_mean = []
    tmrca_err = []

    N_MUS = np.unique(DF.Nmu)
    N_MUS_idxs = np.ones(N_MUS.shape, dtype=bool)

    for idx, N_MU  in enumerate(N_MUS):

        idxs = DF.Nmu == N_MU
        if idxs.sum() == 0:
            N_MUS_idxs[idx] = False
            continue

        dMu = (DF.Sim_Mu[idxs] - DF.Mu[idxs])/DF.Sim_Mu[idxs]
        dMu.sort_values(inplace=True)
        dMu[int(dMu.shape[0]*0.05) : int(dMu.shape[0]*0.95)]

        dTmrca = DF.dTmrca[idxs]/DF.N[idxs]
        dTmrca.sort_values(inplace=True)
        dTmrca = dTmrca[int(dTmrca.shape[0]*0.05) : int(dTmrca.shape[0]*0.95)]

        if mean_or_median == "mean":
            mu_mean.append(np.mean(dMu))
            mu_err.append(np.std(dMu))

            tmrca_mean.append(np.mean(dTmrca))
            tmrca_err.append(np.std(dTmrca))
        else:
            q75, q25 = np.percentile(dMu, [75 ,25])
            mu_err.append((q75 - q25)) #np.std(DF.dTmrca[idxs])
            mu_mean.append(np.median(dMu))
            q75, q25 = np.percentile(dTmrca, [75 ,25])
            tmrca_err.append((q75 - q25)) #np.std(DF.dTmrca[idxs])
            tmrca_mean.append(np.median(dTmrca))

    res = pandas.DataFrame({
        "Nmu" : N_MUS[N_MUS_idxs],
        "dMu_mean" : mu_mean,
        "dMu_err" : mu_err,
        "dTmrca_mean" : tmrca_mean,
        "dTmrca_err" : tmrca_err,
        })
    res = res.sort_values(by='Nmu')
    return res

def plot_data_stat(what, axes, beast=None, tt=None, tt_f=None, lsd=None, lsd_f=None, plot_idxs=None):

    from plot_defaults import shift_point_by_markersize

    axes.grid('on')
    axes.set_xscale('log')

    if what == 'Mu':
        mean = 'dMu_mean'
        err = 'dMu_err'
        title = "Mutation rate deviation"
        ylabel = "Relative mutation rate error, $[\Delta\mu / \mu]$"

    elif what == 'Tmrca':
        mean = 'dTmrca_mean'
        err = 'dTmrca_err'
        title = "Accuracy of Tmrca prediction"
        ylabel = "Relative Tmrca error, $[\Delta\mathrm{T_{mrca}} / \mathrm{N}]$"

    if beast is not None:

        if plot_idxs is None:
            plot_idxs = np.ones(beast.shape[0] ,dtype=bool)

        axes.errorbar(beast["Nmu"].loc[plot_idxs].values, beast[mean].loc[plot_idxs].values, beast[err].loc[plot_idxs].values,
            marker='o',
            markersize=markersize,
            c=beast_col,
            label="Beast")

    if tt is not None:
        x, y = shift_point_by_markersize(axes, tt["Nmu"], tt[mean], -markersize/4)
        if plot_idxs is None:
            plot_idxs = np.ones(x.shape[0] ,dtype=bool)
        axes.errorbar(x[plot_idxs], y[plot_idxs], tt[err]/2,
            fmt='--',
            marker='o',
            markersize=markersize,
            c=tt_col,
            label="TreeTime, original tree")

    if tt_f is not None:
        x, y = shift_point_by_markersize(axes, tt_f["Nmu"], tt_f[mean], +markersize*.75)
        if plot_idxs is None:
            plot_idxs = np.ones(x.shape[0] ,dtype=bool)
        axes.errorbar(x[plot_idxs], y[plot_idxs], (tt_f[err].values/2)[plot_idxs],
            marker='o',
            markersize=markersize,
            markerfacecolor='w',
            markeredgecolor=tt_col,
            mew=1.3,
            c=tt_col, label="TreeTime, FastTree")

    if lsd is not None:
        x, y = shift_point_by_markersize(axes, lsd["Nmu"], lsd[mean], +markersize/2)
        if plot_idxs is None:
            plot_idxs = np.ones(x.shape[0] ,dtype=bool)
        axes.errorbar(x[plot_idxs], y[plot_idxs], (lsd[err].values/2)[plot_idxs],
            fmt='--',
            marker='o',
            markersize=markersize,
            c=lsd_col,
            label="LSd, original tree")

    if lsd_f is not None:
        x, y = shift_point_by_markersize(axes, lsd_f["Nmu"], lsd_f[mean], -markersize*.75)
        if plot_idxs is None:
            plot_idxs = np.ones(x.shape[0] ,dtype=bool)
        axes.errorbar(x[plot_idxs], y[plot_idxs], (lsd_f[err].values/2)[plot_idxs],
            marker='o',
            markersize=markersize,
            markerfacecolor='w',
            markeredgecolor=lsd_col,
            mew=1.3,
            c=lsd_col, label="LSD, FastTree")


    plt.hlines(0, 0, 1)
    axes.legend(loc=0)
    #axes.set_title(title)
    axes.set_ylabel(ylabel, fontsize = label_fs)
    axes.set_xlabel('$\mathrm{N}\cdot\mu$', fontsize = label_fs)
    for label in axes.get_xticklabels():
            label.set_fontsize(tick_fs)
    for label in axes.get_yticklabels():
            label.set_fontsize(tick_fs)


if __name__ == '__main__':

    T_over_N = 10.
    mean_or_median = 'median'
    PLOT_TREETIME = True
    PLOT_LSD = True
    PLOT_BEAST = True
    save_fig = False
    plot_idxs = np.array([1,2,4,6,7,9,10])

    res_dir = "./simulated_data"
    res_lsd = read_lsd_results_dataset('./simulated_data/2017-04-19_lsd_res.csv')
    res_lsd_f = read_lsd_results_dataset('./simulated_data/2017-04-19_lsd_fasttree_res.csv')

    res_tt = read_treetime_results_dataset('./simulated_data/2017-04-19_treetime_res.csv')
    res_tt_f =  read_treetime_results_dataset('./simulated_data/2017-04-19_treetime_fasttree_res.csv')

    beast_logs_dir = os.path.join(res_dir, '2017-04-19_beast')
    beast_trees_dir = os.path.join(res_dir, 'dataset')


    if PLOT_LSD:
        pivot_lsd = None #create_lsd_tt_pivot(res_lsd, T_over_N=T_over_N, mean_or_median=mean_or_median)
        pivot_lsd_f = create_lsd_tt_pivot(res_lsd_f, T_over_N=T_over_N, mean_or_median=mean_or_median)
    else:
        pivot_lsd = None
        pivot_lsd_f = None

    if PLOT_TREETIME:
        pivot_tt = None #create_lsd_tt_pivot(res_tt, T_over_N=T_over_N, mean_or_median=mean_or_median)
        pivot_tt_f =  create_lsd_tt_pivot(res_tt_f, T_over_N=T_over_N, mean_or_median=mean_or_median)
    else:
        pivot_tt = None
        pivot_tt_f = None

    if PLOT_BEAST:
        beast_df = read_all_beast_logs(beast_logs_dir, beast_trees_dir, T_over_N=T_over_N)
        pivot_beast = create_beast_log_pivot(beast_df)
    else:
        pivot_beast=None

    fig = plt.figure(figsize=onecolumn_figsize)
    axes = fig.add_subplot(111)
    plot_data_stat('Mu', axes, beast=pivot_beast, tt=pivot_tt, tt_f=pivot_tt_f, lsd=pivot_lsd, lsd_f=pivot_lsd_f, plot_idxs=plot_idxs)

    if save_fig:
        fig.savefig("./figs/simdata_Mu_TN{}_{}.svg".format(T_over_N, mean_or_median))
        fig.savefig("./figs/simdata_Mu_TN{}_{}.jpg".format(T_over_N, mean_or_median))

    fig = plt.figure(figsize=onecolumn_figsize)
    axes = fig.add_subplot(111)
    plot_data_stat('Tmrca', axes, beast=pivot_beast, tt=pivot_tt, tt_f=pivot_tt_f, lsd=pivot_lsd, lsd_f=pivot_lsd_f, plot_idxs=plot_idxs)

    if save_fig:
        fig.savefig("./figs/simdata_Tmrca_TN{}_{}.svg".format(T_over_N, mean_or_median))
        fig.savefig("./figs/simdata_Tmrca_TN{}_{}.jpg".format(T_over_N, mean_or_median))
