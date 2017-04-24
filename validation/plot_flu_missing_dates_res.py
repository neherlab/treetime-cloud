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

def read_results_dataframe(inf):
    cols = ['File', 'Frac', 'Tmrca', 'Mu', 'Mu_R2', 'Internal_R2', 'Runtime']
    df = pandas.read_csv(inf, names=cols)
    return df

def make_results_pivot(df):

    fs = np.unique(df['Frac'])

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

def plot_results(df, what, label="", axes=None,shift_points=None):
    if what == 'Tmrca':
        ylabel = r'$\mathrm{T}_{mrca}, [year]$'
        title = 'Dependence of $\mathrm{T}_{mrca}$ prediction on the fraction of the known dates'
    elif what == 'Mu':
        ylabel = 'Mutation rate, $[\mathrm{year}^{-1}]$'
        title = 'Dependence of Mutation rate assessment on the fraction of the known dates'
    elif what == 'R2':
        ylabel = 'Molecular clock regression'
        title = 'Quality of the molecular clock regression vs fraction of known dates'

    if axes is None:
        fig = plt.figure(figsize=onecolumn_figsize)
        axes = fig.add_subplot(111)

    if shift_points is not None:
        x, y = shift_point_by_markersize(axes, df.Frac, df[what], shift_points)
    else:
        x, y = df.Frac, df[what]

    axes.ticklabel_format(useOffset=False)
    axes.errorbar(x, y, df[what + '_err'], fmt='o-', label=label,markersize=markersize)
    axes.set_xlabel('Fraction of known dates', fontsize = label_fs)
    axes.set_ylabel(ylabel, fontsize=label_fs)
    axes.set_title(title)
    axes.grid('on')
    axes.legend(loc=0)
    for label in axes.get_xticklabels():
            label.set_fontsize(tick_fs)
    for label in axes.get_yticklabels():
            label.set_fontsize(tick_fs)

def read_dates_stat(inf):
    cols=['tree_name', 'known_frac', 'Tmrca', 'Date_sim','Date_given','dT']
    df = pandas.read_csv(inf, names=cols)
    return df

def make_dates_pivot(df):
    res = pandas.DataFrame()
    fracs = np.unique(df['known_frac'])
    for frac in fracs:
        idxs = df['known_frac'] == frac
        y, x = np.histogram(df['dT'][idxs],bins=20)
        res["F_{}_x".format(frac)] = 0.5*(x[:-1]+x[1:])
        res["F_{}_y".format(frac)] = 1. * y / y.max()
    return res

def plot_date_dists(res, label="", axes=None):

    ## Define FwHM:
    def _fwhm(x, y, points=100):
        interp = interp1d(x, y)
        new_x = np.linspace(interp.x.min(), interp.x.max(), points)
        upper = np.where(interp(np.linspace(interp.x.min(), interp.x.max(), 100)) > 0.5)[0]
        if upper.shape[0]==0:
            return 0
        return new_x[upper[-1]] - new_x[upper[0]]

    if axes is None:
        fig = plt.figure(figsize=onecolumn_figsize)
        axes = fig.add_subplot(111)

    fracs = np.unique([k[:-2] for k in res.columns])

    even = -1
    for frac in fracs:
        even += 1
        if even % 2 != 0:
            continue
        numfrac=float(frac[2:])
        x = frac + "_x"
        y = frac + "_y"
        plt.plot(res[x], res[y], label="{}% Dates known".format(numfrac*100))

    s_x, s_y = [float(k[2:]) for k in fracs], [_fwhm(res[k+"_x"], res[k+"_y"]) for k in fracs]
    #import ipdb; ipdb.set_trace()

    # this is an inset axes over the main axes
    a = plt.axes([.65, .65, .22, .22], axisbg='lightgray',frameon=True)
    a.plot(s_x, s_y, 'o')
    a.set_title('FWHM')
    a.set_xlabel("Fraction of dates known",fontsize=label_fs*.6)
    a.set_ylabel("Leaf date error, $[\mathrm{years}]$",fontsize=label_fs*.6)
    a.get_xaxis().set_ticks([0.1,0.5,0.9])
    a.get_yaxis().set_ticks([0.4,0.7,1.])
    for label in a.get_xticklabels():
            label.set_fontsize(tick_fs*.6)
    for label in a.get_yticklabels():
            label.set_fontsize(tick_fs*.6)

    #plt.xticks([])
    #plt.yticks([])

    axes.set_xlim(-2, 2)
    axes.grid('on')
    axes.legend(loc=2,fontsize=legend_fs)
    axes.set_xlabel(r"Error in $\mathrm{T}_{mrca}$ estimation, $[\mathrm{years}]$", fontsize=label_fs)

    for label in axes.get_xticklabels():
            label.set_fontsize(tick_fs)
    for label in axes.get_yticklabels():
            label.set_fontsize(tick_fs)



if __name__ == '__main__':

    save_fig = True

    work_dir = './flu_H3N2/missing_dates'
    fname_format = 'H3N2_HA_2011_2013_{}seqs_res.csv'
    fname_dates_format = 'H3N2_HA_2011_2013_{}seqs_res.csv_dates.csv'

    df_100 = make_results_pivot(read_results_dataframe(os.path.join(work_dir, fname_format.format(100))))
    df_500 = make_results_pivot(read_results_dataframe(os.path.join(work_dir, fname_format.format(500))))

    dates_100 = make_dates_pivot(read_dates_stat( os.path.join(work_dir, fname_dates_format.format(100))))

    fig = plt.figure(figsize=onecolumn_figsize)
    ax_dates = fig.add_subplot(111)
    plot_date_dists(dates_100,axes=ax_dates)
    if save_fig:
        fig.savefig("./figs/fluH3N2_missingDates_leafDates.svg")
        fig.savefig("./figs/fluH3N2_missingDates_leafDates.jpg")

    fig = plt.figure(figsize=onecolumn_figsize)
    ax_tmrca = fig.add_subplot(111)
    plot_results(df_100, 'Tmrca', label='100 nodes tree', axes=ax_tmrca, shift_points=+markersize/2)
    plot_results(df_500, 'Tmrca', label='500 nodes tree', axes=ax_tmrca, shift_points=-markersize/2)
    if save_fig:
        fig.savefig("./figs/fluH3N2_missingDates_Tmrca.svg")
        fig.savefig("./figs/fluH3N2_missingDates_Tmrca.jpg")


    fig = plt.figure(figsize=onecolumn_figsize)
    ax_mu = fig.add_subplot(111)
    plot_results(df_100, 'Mu', label='100 nodes tree', axes=ax_mu, shift_points=+markersize/2)
    plot_results(df_500, 'Mu', label='500 nodes tree', axes=ax_mu, shift_points=-markersize/2)
    if save_fig:
        fig.savefig("./figs/fluH3N2_missingDates_Mu.svg")
        fig.savefig("./figs/fluH3N2_missingDates_Mu.jpg")

    fig = plt.figure(figsize=onecolumn_figsize)
    ax_mu = fig.add_subplot(111)
    plot_results(df_100, 'R2', label='100 nodes tree', axes=ax_mu, shift_points=+markersize/2)
    plot_results(df_500, 'R2', label='500 nodes tree', axes=ax_mu, shift_points=-markersize/2)
    if save_fig:
        fig.savefig("./figs/fluH3N2_missingDates_MolClockR2.svg")
        fig.savefig("./figs/fluH3N2_missingDates_MolClockR2.jpg")
