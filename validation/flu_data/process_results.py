import pandas
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import os, sys
import shutil
plt.ion()
from Bio import Phylo
from analysis import dates_from_flu_tree

def read_lsd_dataset(fname):

    lsd_cols = ['File', 'N', 'Tmrca_sim', 'mu_sim', 'Runtime', 'objective']
    lsd_df = pandas.read_csv(fname, names=lsd_cols)
    return lsd_df

def read_treetime_dataset(fname):

    cols = ['File', 'N', "Tmrca_sim", "mu_sim", "R2_leaves", "R2_internal", "Runtime"]
    df = pandas.read_csv(fname, names=cols)

    return df

def read_beast_log(logfile, nearest_leaf_date):
    df = pandas.read_csv(logfile, delimiter='\t', skiprows=3)
    df['Tmrca_stat']  = nearest_leaf_date  - df['Tmrca_stat']
    return df.iloc[-500:]


def beast_logs_stat(logsdir='./beast_out_cluster_1', treedir='./subtrees'):

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
        dates = dates_from_flu_tree(treepath)

        df = read_beast_log(logpath, np.max(dates.values()))
        if df is None or df.shape[0] < 300:
            continue
        new_Ns = Phylo.read(treepath, 'newick').count_terminals()
        if new_Ns > 1e3:
            continue
        Ns.append(new_Ns)

        LH_mean.append(df['likelihood'].mean())
        LH_err.append(df['likelihood'].std())

        Mu_mean.append(df['meanRate'].mean())
        Mu_err.append(df['meanRate'].std()) # * Mu_mean[-1] / 2.)

        Tmrca_mean.append(df['Tmrca_stat'].mean())
        Tmrca_err.append(df['Tmrca_stat'].std())

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

def plot_res(what, tt=None, lsd=None, beast=None, save=True, suffix=None, **kwargs):


    if what == 'Tmrca':
        mean = 'Tmrca_mean'
        err = 'Tmrca_err'
        title = "Estimated Tmrca as function of sample size\nLSD params: -{}".format(suffix)
        ylabel = "Tmrca"

    elif what == "Mu":
        mean = 'Mu_mean'
        err =  'Mu_err'
        title =  "Estimated Mutation rate as function of sample size\nLSD params: -{}".format(suffix)
        ylabel = "Mutation rate, #/year"

    fig = plt.figure()
    axes = fig.add_subplot(111)

    if tt is not None:
        axes.errorbar(tt['Ns'], tt[mean], tt[err], markersize=10, marker='o', c='b', label='TreeTime')

    if lsd is not None:
        axes.errorbar(lsd['Ns'], lsd[mean], lsd[err], markersize=10, marker='o', c='g', label='LSD')

    if beast is not None:
        axes.errorbar(beast['Ns'], beast[mean], beast[err], markersize=10, marker='o', c='r', label='BEAST')

    axes.grid()
    axes.legend(loc=0)
    axes.set_xscale('log')

    axes.set_ylabel(ylabel)
    axes.set_xlabel("Tree size, #sequences")
    axes.set_title(title)

    if save:
        fig.savefig("./figs/{}_LSD_{}.svg".format(ylabel, suffix))
        fig.savefig("./figs/{}_LSD_{}.jpg".format(ylabel, suffix))


if __name__ == "__main__":

    #fig = plt.figure()


    #plot_beast_logs(axes, what='meanRate')

#    suffix = sys.argv[1]
#
#
    tt_df = treetime_dataset_stat(read_treetime_dataset('./H3N2_HA_1980_2015_RES_treetime.csv'))
    lsd_df = lsd_dataset_stat(read_lsd_dataset('./H3N2_HA_1980_2015_RES_lsd.csv'))
    beast = beast_logs_stat(logsdir='./beast_out_cluster_2')

    plot_res('Tmrca', tt_df, lsd_df, beast, save=True)

    plot_res('Mu', tt_df, lsd_df, beast, save=True)
#
#    Ns = np.unique(np.concatenate((tt_df["N"], lsd_df['N'])))
#
#    Tmrca_mean_tt = []
#    Tmrca_err_tt = []
#    Mu_mean_tt = []
#    Mu_err_tt = []
#
#    Tmrca_mean_lsd = []
#    Tmrca_err_lsd = []
#    Mu_mean_lsd = []
#    Mu_err_lsd = []
#
#    Runtime_mean_lsd = []
#    Runtime_mean_tt = []
#    Runtime_err_lsd = []
#    Runtime_err_tt = []
#
#
#    for N in Ns:
#
#        Nidx_tt = tt_df["N"] == N
#        Tmrca_mean_tt.append((N, tt_df[Nidx_tt]["Tmrca_sim"].mean()))
#        Tmrca_err_tt.append((N, tt_df[Nidx_tt]["Tmrca_sim"].std()))
#        Mu_mean_tt.append((N, tt_df[Nidx_tt]["mu_sim"].mean()))
#        Mu_err_tt.append((N, tt_df[Nidx_tt]["mu_sim"].std()))
#        Runtime_mean_tt.append((N, tt_df[Nidx_tt]["Runtime"].mean()))
#        Runtime_err_tt .append((N, tt_df[Nidx_tt]["Runtime"].std()))
#
#        Nidx_lsd = lsd_df["N"] == N
#        Tmrca_mean_lsd.append((N, lsd_df[Nidx_lsd]["Tmrca_sim"].mean()))
#        Tmrca_err_lsd.append((N, lsd_df[Nidx_lsd]["Tmrca_sim"].std()))
#        Mu_mean_lsd.append((N, lsd_df[Nidx_lsd]["mu_sim"].mean()))
#        Mu_err_lsd.append((N, lsd_df[Nidx_lsd]["mu_sim"].std()))
#        Runtime_mean_lsd.append((N, lsd_df[Nidx_lsd]["Runtime"].mean()))
#        Runtime_err_lsd .append((N, lsd_df[Nidx_lsd]["Runtime"].std()))
#
#
#    Tmrca_mean_tt = np.array(Tmrca_mean_tt)
#    Tmrca_err_tt = np.array(Tmrca_err_tt)
#    Mu_mean_tt = np.array(Mu_mean_tt)
#    Mu_err_tt = np.array (Mu_err_tt)
#
#    Tmrca_mean_lsd = np.array(Tmrca_mean_lsd)
#    Tmrca_err_lsd = np.array(Tmrca_err_lsd)
#    Mu_mean_lsd = np.array(Mu_mean_lsd)
#    Mu_err_lsd =  np.array (Mu_err_lsd)
#
#    Runtime_mean_lsd = np.array(Runtime_mean_lsd)
#    Runtime_mean_tt  = np.array(Runtime_mean_tt)
#    Runtime_err_lsd  = np.array(Runtime_err_lsd)
#    Runtime_err_tt   = np.array(Runtime_err_tt)
#
#
#    plt.figure()
#    plt.errorbar(Tmrca_err_tt[:, 0], Tmrca_mean_tt[:, 1], Tmrca_err_tt[:, 1], fmt='o-', label='TreeTime')
#    plt.errorbar(Tmrca_err_lsd[:, 0], Tmrca_mean_lsd[:, 1], Tmrca_err_lsd[:, 1], fmt='o-', label='LSD')
#    plt.grid()
#    plt.legend()
#    plt.ylabel("Tmrca")
#    plt.title("Estimated Tmrca as function of sample size\nLSD params: -{}".format(suffix))
#    plt.xlabel("Tree size, #sequences")
#    plt.xscale("log")
#    plt.savefig("./figs/Tmrca_LSD_{}.svg".format(suffix))
#    plt.savefig("./figs/Tmrca_LSD_{}.jpg".format(suffix))
#
#
#    plt.figure()
#    plt.errorbar(Mu_err_tt[:, 0],  Mu_mean_tt[:, 1],  Mu_err_tt[:, 1], fmt='o-', label='TreeTime')
#    plt.errorbar(Mu_err_lsd[:, 0], Mu_mean_lsd[:, 1], Mu_err_lsd[:, 1], fmt='o-', label='LSD')
#    plt.grid()
#    plt.legend()
#    plt.ylabel("Mutation rate")
#    plt.title("Estimated Mutation rate as function of sample size\nLSD params: -{}".format(suffix))
#    plt.xlabel("Tree size, #sequences")
#    plt.xscale("log")
#    plt.savefig("./figs/Mu_LSD_{}.svg".format(suffix))
#    plt.savefig("./figs/Mu_LSD_{}.jpg".format(suffix))
#
#
#
#    plt.figure()
#    plt.errorbar(Runtime_err_tt[:, 0],  Runtime_mean_tt[:, 1],  Runtime_err_tt[:, 1], fmt='o-', label='TreeTime')
#    plt.errorbar(Runtime_err_lsd[:, 0], Runtime_mean_lsd[:, 1], Runtime_err_lsd[:, 1], fmt='o-', label='LSD')
#    plt.grid()
#    plt.legend()
#    plt.ylabel("Runtime, arb.u.")
#    plt.title("Estimated Runtime as function of sample size")
#    plt.xlabel("Sample size, MAX # sequences per year")
#
#
#    ft  = "./LSD/H3N2_HA_1980_2015_NA_20_23.txt.newick"
#    shutil.copyfile(ft, "./figs/tree_LSD_{}.nwk".format(suffix))



