import numpy as np
import datetime
import os, copy
import matplotlib.pyplot as plt
import pandas
from scipy.interpolate import interp1d
from plot_defaults import *

def read_dates_stat(inf):
    cols=['tree_name', 'known_frac', 'Tmrca', 'Date_sim','Date_given','dT']
    df = pandas.read_csv(inf, names=cols,header=0)
    return df

def make_dates_pivot(df):
    res = pandas.DataFrame()
    fracs = np.unique(df['known_frac'])
    print (fracs)
    for frac in fracs:
        idxs = df['known_frac'] == frac
        y, x = np.histogram(df['dT'][idxs],bins=20)
        res["F_{}_x".format(frac)] = 0.5*(x[:-1]+x[1:])
        res["F_{}_y".format(frac)] = 1. * y / y.max()
    return res

def plot_date_dists(res, label="", figname=None, idx=[0.05, 0.2, 0.4, 0.8]):

    ## Define FwHM:
    def _fwhm(x, y, points=100):
        interp = interp1d(x, y)
        new_x = np.linspace(interp.x.min(), interp.x.max(), points)
        upper = np.where(interp(np.linspace(interp.x.min(), interp.x.max(), 100)) > 0.5)[0]
        if upper.shape[0]==0:
            return 0
        return new_x[upper[-1]] - new_x[upper[0]]

    fig = plt.figure(figsize=onecolumn_figsize)
    axes = fig.add_subplot(111)

    fracs = np.unique([k[:-2] for k in res.columns])

    even = -1
    for frac in fracs:
        numfrac=float(frac[2:])
        if numfrac not in idx:
            continue
        x = frac + "_x"
        y = frac + "_y"
        axes.plot(res[x], res[y], label="{}% dates known".format(numfrac*100), lw=2)

    s_x, s_y = [float(k[2:]) for k in fracs], [_fwhm(res[k+"_x"], res[k+"_y"]) for k in fracs]

    # this is an inset axes over the main axes
    a = plt.axes([.61, .6, .27, .27], axisbg='lightgray',frameon=True)
    a.plot(s_x, s_y, 'o', markersize=markersize*0.8)
    a.set_title('FWHM')
    a.set_xlabel("Fraction of dates known",fontsize=label_fs*.8)
    a.set_ylabel("Leaf date error, $[\mathrm{Years}]$",fontsize=label_fs*.8)
    a.get_xaxis().set_ticks([0.0,0.25,0.5,0.75,1.])
    a.get_yaxis().set_ticks([0.4,0.6,0.8,1.])
    for label in a.get_xticklabels():
            label.set_fontsize(tick_fs*.8)
    for label in a.get_yticklabels():
            label.set_fontsize(tick_fs*.8)

    #plt.xticks([])
    #plt.yticks([])

    axes.set_xlim(-2, 2)
    axes.set_ylim(0, 1.4)
    axes.grid('on')
    axes.legend(loc=2,fontsize=legend_fs)
    axes.set_xlabel(r"Error in leaf dates reconstruction, $[\mathrm{Years}]$", fontsize=label_fs)

    for label in axes.get_xticklabels():
            label.set_fontsize(tick_fs)
    for label in axes.get_yticklabels():
            label.set_fontsize(tick_fs)

    if figname is not None:
        for fmt in formats:
            fig.savefig("{}.{}".format(figname, fmt))

if __name__ == '__main__':

    SAVE_FIG = False

    ##
    ## file with the reconstruction results:
    ##
    leaf_dates_file = './flu_H3N2/missing_dates/treetime_dates_res.csv'

    ##
    ## Read the leaf dates file, make pivot table:
    ##
    print ("Plotting the results")
    dataframe = make_dates_pivot(read_dates_stat(leaf_dates_file))

    ##
    ## Plot the leaf dates distribution
    ##
    plot_date_dists(dataframe,
        figname="./figs/fluH3N2_missingDates_leafDates" if SAVE_FIG else None)

