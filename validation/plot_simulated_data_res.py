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
    df['dTmrca'] = -(df['Sim_Tmrca'] - df['Tmrca'])
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
    df['dTmrca'] = -(df['Sim_Tmrca'] - df['Tmrca'])
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

#        import ipdb; ipdb.set_trace()
        dMu = -(DF.Sim_mu[idxs] - DF.mu[idxs])/DF.Sim_mu[idxs]
        dMu.sort_values(inplace=True)
        #dMu = dMu[int(dMu.shape[0]*0.05) : int(dMu.shape[0]*0.95)]

        dTmrca = DF.dTmrca[idxs]/DF.N[idxs]
        dTmrca.sort_values(inplace=True)
        #dTmrca = dTmrca[int(dTmrca.shape[0]*0.05) : int(dTmrca.shape[0]*0.95)]

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
    fig = plt.figure(figsize=onecolumn_figsize)
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(NUM_COLORS)])

    for MU in N_MUS:
        #ax.plot(df["T"]/df['N'], df.dTmrca / (2016.5 - df.Tmrca), 'o', label="MU = " + str(MU), alpha=0.7)
        ax.plot(df["T"]/df['N']*(1+np.random.normal(size=len(df))*0.1), df.dTmrca / df["N"], 'o', label="MU = " + str(MU), alpha=0.6)

    ax.legend(loc=0, fontsize=legend_fs)
    plt.xlabel("Total evolution time,  $T/N$")
    plt.ylabel("Relative error in $T_{mrca}$ estimation, $\Delta T_{mrca} / N$")

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

    fig = plt.figure(figsize=onecolumn_figsize)
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
    plt.legend(loc=0, fontsize=legend_fs)

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

        dTmrca.append(-(Sim_Tmrca[-1] - Tmrca[-1]))
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

        dMu = -(DF.Sim_Mu[idxs] - DF.Mu[idxs])/DF.Sim_Mu[idxs]
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
        title = "Clock rate deviation"
        ylabel = "Relative clock rate error, $[\Delta\mu / \mu]$"

    elif what == 'Tmrca':
        mean = 'dTmrca_mean'
        err = 'dTmrca_err'
        title = "Accuracy of Tmrca prediction"
        ylabel = "Relative $T_{mrca}$ error, $[\Delta\mathrm{T_{mrca}} / \mathrm{N}]$"

    if beast is not None:
        if plot_idxs is None:
            beast_plot_idxs = np.ones(beast.shape[0] ,dtype=bool)
        else:
            beast_plot_idxs = plot_idxs

        axes.errorbar(beast["Nmu"].loc[beast_plot_idxs].values, beast[mean].loc[beast_plot_idxs].values, beast[err].loc[beast_plot_idxs].values,
            marker='o',
            markersize=markersize,
            c=beast_color,
            label="Beast")

    if tt is not None:
        x, y = shift_point_by_markersize(axes, tt["Nmu"], tt[mean], +markersize*.75)
        if plot_idxs is None:
            tt_plot_idxs = np.ones(x.shape[0] ,dtype=bool)
        else:
            tt_plot_idxs = plot_idxs

        axes.errorbar(x[tt_plot_idxs], y[tt_plot_idxs], (tt[err].values/2)[tt_plot_idxs],
            fmt='-',
            marker='o',
            markersize=markersize,
            #markerfacecolor='w',
            markeredgecolor=tt_color,
            mew=1.3,
            c=tt_color, label="TreeTime, original tree")


    if tt_f is not None:
        x, y = shift_point_by_markersize(axes, tt_f["Nmu"], tt_f[mean], +markersize*.75)
        if plot_idxs is None:
            ttf_plot_idxs = np.ones(x.shape[0] ,dtype=bool)
        else:
            ttf_plot_idxs = plot_idxs

        axes.errorbar(x[ttf_plot_idxs], y[ttf_plot_idxs], (tt_f[err].values/2)[ttf_plot_idxs],
            marker='o',
            markersize=markersize,
            markerfacecolor='w',
            markeredgecolor=tt_color,
            mew=1.3,
            c=tt_color, label="TreeTime, FastTree")

    if lsd is not None:
        x, y = shift_point_by_markersize(axes, lsd["Nmu"], lsd[mean], +markersize/2)
        if plot_idxs is None:
            lsd_plot_idxs = np.ones(x.shape[0] ,dtype=bool)
        else:
            lsd_plot_idxs = plot_idxs

        axes.errorbar(x[lsd_plot_idxs], y[lsd_plot_idxs], (lsd[err].values/2)[lsd_plot_idxs],
            fmt='--',
            marker='o',
            markersize=markersize,
            c=lsd_color,
            label="LSd, original tree")

    if lsd_f is not None:
        x, y = shift_point_by_markersize(axes, lsd_f["Nmu"], lsd_f[mean], -markersize*.75)
        if plot_idxs is None:
            lsdf_plot_idxs = np.ones(x.shape[0] ,dtype=bool)
        else:
            lsdf_plot_idxs = plot_idxs

        axes.errorbar(x[lsdf_plot_idxs], y[lsdf_plot_idxs], (lsd_f[err].values/2)[lsdf_plot_idxs],
            marker='o',
            markersize=markersize,
            markerfacecolor='w',
            markeredgecolor=lsd_color,
            mew=1.3,
            c=lsd_color, label="LSD, FastTree")


    plt.hlines(0, 0, 1)
    axes.legend(loc=1,fontsize=legend_fs)
    #axes.set_title(title)
    axes.set_ylabel(ylabel, fontsize = label_fs)
    axes.set_xlabel('Diversity, $\mathrm{N}\cdot\mu$', fontsize = label_fs)
    for label in axes.get_xticklabels():
            label.set_fontsize(tick_fs)
    for label in axes.get_yticklabels():
            label.set_fontsize(tick_fs)

def get_beast_tree_from_file(beast_file):
    import dendropy
    import StringIO
    trees = dendropy.TreeList.get(path=beast_file, schema="nexus")
    tree = trees[-1]
    out = StringIO.StringIO()
    tree.write_to_stream(out, 'newick')
    biotree = Phylo.read(StringIO.StringIO(out.getvalue()), 'newick')
    return biotree


def corr_points(basename, beast_dir=None):

    def initialize_splits(tree):

        all_tips = set(tree.get_terminals())
        for c in tree.find_clades(order='postorder'):
            if c.is_terminal():
                c.tips = [c.name]
            else:
                c.tips = set(sum([list(k.tips) for k in c.clades],[]))
            c.non_tips = set.difference(all_tips, c.tips)

        return tree

    def get_beast_tree_from_file(beast_file):
        print ("reading BEAST tree from file: " + beast_file)
        if beast_file is None:
            return None
        import dendropy
        import StringIO
        trees = dendropy.TreeList.get_from_path(beast_file, schema="nexus")
        tree = trees[-1]
        out = StringIO.StringIO()
        tree.write_to_stream(out, 'newick')
        biotree = Phylo.read(StringIO.StringIO(out.getvalue()), 'newick')
        return biotree

    def basename_to_beast_file(basename, beast_dir):
        if beast_dir is None:
            return None
        return os.path.join(beast_dir, os.path.split(basename)[-1] + '.trees.txt')


    treetime_tree = initialize_splits(Phylo.read(basename + '.treetrime.ft.nwk', 'newick'))
    original_tree = initialize_splits(Phylo.read(basename + '.nwk', 'newick'))
    fasttree_tree = initialize_splits(Phylo.read(basename + '.ft.nwk', 'newick'))
    beast_tree = get_beast_tree_from_file(basename_to_beast_file(basename, beast_dir))
    if beast_tree is not None:
        beast_tree = initialize_splits(beast_tree)

    tt_corr = []
    ft_corr = []
    bt_corr = []

    for orig_split in original_tree.find_clades():
        for c1 in treetime_tree.find_clades():
            if c1.tips == orig_split.tips:
                tt_corr.append((orig_split.branch_length, c1.branch_length))
                break
        for c2 in fasttree_tree.find_clades():
            if c2.tips == orig_split.tips:
                ft_corr.append((orig_split.branch_length, c2.branch_length))
        if beast_tree is not None:
            for c3 in beast_tree.find_clades():
                if c3.tips == orig_split.tips:
                    bt_corr.append((orig_split.branch_length, c3.branch_length))

    return tt_corr, ft_corr, bt_corr

def correlation_dataset(root_dir, beast_root_dir, Mu=['0.0001'], Ts=['50'], **kwargs):

    basenames = [os.path.join(root_dir, k[:-17]) for k in os.listdir(root_dir)
        if 'treetrime.ft.nwk' in k
        and 'Ts{}'.format(Ts[0]) in k
        and 'Mu{}'.format(Mu[0]) in k]

    #basenames = [basenames[0]]

    cors = [corr_points(basename, beast_root_dir) for basename in basenames]
    tt_corr = np.array(sum([k[0] for k in cors], []))
    ft_corr = np.array(sum([k[1] for k in cors], []))
    bt_corr = np.array(sum([k[2] for k in cors], []))
    tt_corr[:, 1] /= float(Mu[0])
    ft_corr[:, 1] /= float(Mu[0])
    return tt_corr, ft_corr, bt_corr

def plot_correlation(tt_corr, ft_corr, bt_corr, axes=None, include_fast_tree=True, **kwargs):

    tt_corrcoeff = np.corrcoef(tt_corr[:, 0], tt_corr[:, 1])
    bt_corrcoeff = np.corrcoef(bt_corr[:, 0], bt_corr[:, 1])
    ft_corrcoeff = np.corrcoef(ft_corr[:, 0], ft_corr[:, 1])

    if axes is None:
        fig = plt.figure(figsize=onecolumn_figsize)
        axes = fig.add_subplot(111)

    if include_fast_tree:
        axes.plot(ft_corr[:, 0], ft_corr[:, 1], 'o',
                alpha=0.5,
                marker='o',
                markersize=markersize,
                c=ft_color,
                label="ML reconstruction (FastTree). Correlation coefficient: " + format(ft_corrcoeff[0, 1], '.3f'))

    axes.plot(bt_corr[:, 0], bt_corr[:, 1], 'o',
            alpha=0.5,
            marker='o',
            markersize=markersize,
            c=beast_color,
            label="BEAST estimation. Correlation coefficient: "+ format(bt_corrcoeff[0, 1], '.3f'))

    axes.plot(tt_corr[:, 0], tt_corr[:, 1], 'o',
            alpha=0.5,
            marker='o',
            markersize=markersize,
            c=tt_color,
            label="TreeTime estimation. Correlation coefficient: " + format(tt_corrcoeff[0, 1], '.3f'))

    axes.legend(loc=0,fontsize=legend_fs)
    axes.grid('on')
    #axes.set_title(title)
    axes.set_ylabel('Reconstructed branch length, [$\mathrm{Years}$]', fontsize = label_fs)
    axes.set_xlabel('Simulated branch length, [$\mathrm{Years}$]', fontsize = label_fs)
    for label in axes.get_xticklabels():
            label.set_fontsize(tick_fs)
    for label in axes.get_yticklabels():
            label.set_fontsize(tick_fs)
    axes.set_xlim(0.1, 200)
    axes.set_ylim(0.1, 200)
    #axes.set_xscale('log')
    #axes.set_yscale('log')

if __name__ == '__main__':

    T_over_N = 2.
    mean_or_median = 'median'

    PLOT_SIM_RESULTS = True
    PLOT_CORRELATION = False

    PLOT_TREETIME = True
    PLOT_LSD = False
    PLOT_BEAST = False
    save_fig = False
    plot_idxs = np.array([1,2,4,6,7,9,10])

    if PLOT_SIM_RESULTS:
        res_dir = "./simulated_data"
        res_lsd = read_lsd_results_dataset('./simulated_data/2017-05-31_lsd_res.csv')
        res_lsd_f = read_lsd_results_dataset('./simulated_data/2017-05-31_lsd_fasttree_res.csv')

        res_tt = read_treetime_results_dataset('./simulated_data/2017-06-11_treetime_fasttree_res_use_input_branch_false.csv')
        res_tt_f =  read_treetime_results_dataset('./simulated_data/2017-06-11_treetime_fasttree_res.csv')

        beast_logs_dir = os.path.join(res_dir, '2017-05-16_beast')
        beast_trees_dir = os.path.join(res_dir, 'dataset')


        if PLOT_LSD:
            pivot_lsd = None #create_lsd_tt_pivot(res_lsd, T_over_N=T_over_N, mean_or_median=mean_or_median)
            pivot_lsd_f = create_lsd_tt_pivot(res_lsd_f, T_over_N=T_over_N, mean_or_median=mean_or_median)
        else:
            pivot_lsd = None
            pivot_lsd_f = None

        if PLOT_TREETIME:
            pivot_tt = create_lsd_tt_pivot(res_tt, T_over_N=T_over_N, mean_or_median=mean_or_median)
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
        fig.text(0.15, 0.85, '$\mathrm{\mu}$ overestimated', fontsize=tick_fs)
        fig.text(0.15, 0.15, '$\mathrm{\mu}$ underestimated', fontsize=tick_fs)

        if save_fig:
            fig.savefig("./figs/simdata_Mu_TN{}_{}.svg".format(T_over_N, mean_or_median))
            fig.savefig("./figs/simdata_Mu_TN{}_{}.png".format(T_over_N, mean_or_median))
            fig.savefig("./figs/simdata_Mu_TN{}_{}.pdf".format(T_over_N, mean_or_median))

        fig = plt.figure(figsize=onecolumn_figsize)
        axes = fig.add_subplot(111)
        plot_data_stat('Tmrca', axes, beast=pivot_beast, tt=pivot_tt, tt_f=pivot_tt_f, lsd=pivot_lsd, lsd_f=pivot_lsd_f, plot_idxs=plot_idxs)
        fig.text(0.15, 0.85, '$\mathrm{T_{mrca}}$ too late', fontsize=tick_fs)
        fig.text(0.15, 0.15, '$\mathrm{T_{mrca}}$ too early', fontsize=tick_fs)

        if save_fig:
            fig.savefig("./figs/simdata_Tmrca_TN{}_{}.svg".format(T_over_N, mean_or_median))
            fig.savefig("./figs/simdata_Tmrca_TN{}_{}.png".format(T_over_N, mean_or_median))
            fig.savefig("./figs/simdata_Tmrca_TN{}_{}.pdf".format(T_over_N, mean_or_median))


    if PLOT_CORRELATION:
        INCLUDE_FAST_TREE=True
        root_dir = beast_trees_dir
        beast_root_dir = beast_logs_dir
        tt_corr, ft_corr, bt_corr = correlation_dataset(root_dir, beast_root_dir)
        fig = plt.figure(figsize=onecolumn_figsize)
        axes = fig.add_subplot(111)
        plot_correlation(tt_corr, ft_corr, bt_corr, axes, include_fast_tree=INCLUDE_FAST_TREE)
        if save_fig:
            ft = '_ft' if INCLUDE_FAST_TREE else ""
            fig.savefig("./figs/simdata_BranchLenCorr{}.svg".format(ft))
            fig.savefig("./figs/simdata_BranchLenCorr{}.png".format(ft))
            fig.savefig("./figs/simdata_BranchLenCorr{}.pdf".format(ft))
