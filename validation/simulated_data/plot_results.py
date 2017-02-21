import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors

import pandas

def read_treetime_results_file(fname):
    """
    Read results of the TreeTime simulations

    Args:
     - fname: path to the input file

    Returns:
     - df: Table of results as pandas data-frame
    """

    columns = ['File', 'Sim_Tmrca', 'Tmrca', 'mu', 'R', 'R2_int']

    df = pandas.read_csv(fname, names=columns)
    df = df[len(df["File"]) > 10]

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

def read_lsd_results_file(fname):
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

def plot_raw_data(df):
    """
    Plot the raw simulations data, clustered by the $N\mu$ value
    """

    # number of mutation classes
    MUS = np.unique(df.Nmu)

    #  set colors
    NUM_COLORS = len(MUS)
    cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)

    # plot results
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(NUM_COLORS)])

    for MU in MUS:
        #ax.plot(df["T"]/df['N'], df.dTmrca / (2016.5 - df.Tmrca), 'o', label="MU = " + str(MU), alpha=0.7)
        ax.plot(df["T"]/df['N']*(1+np.random.normal(size=len(df))*0.1), df.dTmrca / df["N"], 'o', label="MU = " + str(MU), alpha=0.6)

    ax.legend()
    plt.xlabel("Total evolution time,  $T/N$")
    plt.ylabel("Relative error in Tmrca estimation, $\Delta T_{mrca} / N$")


def plot_mutation_rate_distributions(df, TN_min=3, TN_max=10, plot_title=""):
    """
    TODO
    """

    # number of mutation classes
    MUS = np.unique(df.Nmu)


    #  set colors
    NUM_COLORS = len(MUS)
    cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)

    DF = df[ (df['T']/df['N'] > TN_min) & (df['T']/df['N'] < TN_max) ]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(NUM_COLORS)])




    for idx,MU  in enumerate(MUS):
        idxs = DF.Nmu == MU
        print (MU, np.sum(idxs))
        if idxs.sum() == 0:
            continue

        # plot histogram
        dMu = (DF.Sim_mu[idxs] - DF.mu[idxs])/DF.Sim_mu[idxs]
        plt.hist(dMu, label="$N\cdot\mu = $" + str(MU), alpha=.5, bins=20)

    ax.set_title(plot_title)
    plt.legend()

def plot_mutation_rate_comparison(treetime, lsd,
                        treetime_fasttree=None, lsd_fasttree=None,
                        TN_min=3, TN_max=10,
                        mu_or_Tmrca='mu',
                        mean_or_median="mean",
                        plot_title="",
                        label="original tree"):
    """
    TODO
    """

    def plot_mu(ax, df, color, linetype, label):
        # filter out Nmu values

        markeredgewidth=1
        markeredgecolor=color
        markerfacecolor='None'
        lw=2

        DF = df[ (df['T']/df['N'] > TN_min) & (df['T']/df['N'] < TN_max) ]

        # number of mutation classes
        MUS = np.unique(df.Nmu)

        mu_stds = np.zeros((len(MUS), 2))
        mu_means = np.zeros((len(MUS), 2))

        for idx,MU  in enumerate(MUS):
            idxs = DF.Nmu == MU
            print (MU, np.sum(idxs))
            if idxs.sum() == 0:
                continue

            if mu_or_Tmrca == 'mu':
                dMu = (DF.Sim_mu[idxs] - DF.mu[idxs])/DF.Sim_mu[idxs]
            else:
                dMu = DF.dTmrca[idxs]/DF.N[idxs]

            mu_stds[idx, 0] = MU
            mu_means[idx,0] = MU
            if mean_or_median == "mean":
                mu_means[idx, 1] = np.mean(dMu)
                mu_stds[idx, 1] = np.std(dMu)
            else:
                q75, q25 = np.percentile(dMu, [75 ,25])
                mu_stds[idx, 1] =  (q75 - q25) #np.std(DF.dTmrca[idxs])
                mu_means[idx,1] = np.median(dMu)


        ax.plot(mu_means[:, 0], mu_means[:, 1], 'o'+linetype, markersize=12,
                                                c=color,
                                                markeredgewidth=markeredgewidth,
                                                markerfacecolor=markerfacecolor,
                                                markeredgecolor=markeredgecolor,
                                                label=label)
        ax.errorbar(mu_means[:, 0], mu_means[:, 1], yerr=mu_stds[:,1]/2., lw=lw/2, fmt=linetype, c=color)
        ax.hlines(0, mu_means[:, 0].min(), mu_means[:, 0].max())


    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_mu(ax, treetime, 'b', '-', "TreeTime, Original tree")
    plot_mu(ax, lsd,      'r', '-', "LSD, Original tree")
    if treetime_fasttree is not None:
        plot_mu(ax, treetime_fasttree, 'b', '--', "TreeTime, Reconstructed tree")
    if lsd_fasttree is not None:
        plot_mu(ax, lsd_fasttree,      'r', '--', "LSD, Reconstructed tree")

    plt.xscale('log')
    ax.set_title(plot_title)
    plt.xlabel('N$\cdot\mu$')
    if mu_or_Tmrca == 'mu':
        plt.ylabel('Relative error $\Delta\mu / \mu$')
    else:
        plt.ylabel('Relative error $\Delta Tmrca / N$')

    plt.legend()
    plt.grid()


if __name__ == '__main__':

    res_lsd = read_lsd_results_file('./accuracy_5__lsd_res.csv')
    res_tt = read_treetime_results_file('./accuracy_5__treetime_res.csv')
    res_lsd_f = read_lsd_results_file('./accuracy_5__lsd_fasttree_res.csv')
    res_tt_f = read_treetime_results_file('./accuracy_5__treetime_fasttree_res.csv')

    plot_raw_data(res_tt)
    plot_raw_data(res_lsd)

    plot_mutation_rate_distributions(res_tt, plot_title="Mutation rate distribution, TreeTime")
    plot_mutation_rate_distributions(res_lsd, plot_title="Mutation rate distribution, LSD")

    plot_mutation_rate_comparison(res_tt, res_lsd, res_tt_f, res_lsd_f,
        mu_or_Tmrca='Tmrca', mean_or_median='median', plot_title="Tmrca inferrence by TreeTime and LSD")

    plot_mutation_rate_comparison(res_tt, res_lsd, res_tt_f, res_lsd_f,
        mu_or_Tmrca='mu', mean_or_median='median', plot_title="Mutation rate inferrence by TreeTime and LSD")


