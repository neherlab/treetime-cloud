"""
Treetime validation on the simulated dataset.
Plot the accuracy of the internal branch lengths reconstruction.

Requires: all trees simulated by FFpopSim, reconstructed with FastTree and TreeTime
should be stored in one directory. Trees reconstructed by BEAST should be also available
and stored (optionally) in another directory. Normally, this is provided by
runing the 'generate_simulated_data_submit.py' script with relevant parameters

The script will scan all trees, identify the similar branches and determine the
accuracy of the branch length reconstruction.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import os, sys
import pandas
from Bio import Phylo
import dendropy
import StringIO

import utility_functions_beast as beast_utils
import utility_functions_simulated_data as sim_utils

from plot_defaults import *


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

def plot_correlation(tt_corr, ft_corr, bt_corr, include_fast_tree=True, figname=None, **kwargs):

    tt_corrcoeff = np.corrcoef(tt_corr[:, 0], tt_corr[:, 1])
    bt_corrcoeff = np.corrcoef(bt_corr[:, 0], bt_corr[:, 1])
    ft_corrcoeff = np.corrcoef(ft_corr[:, 0], ft_corr[:, 1])

    fig = plt.figure(figsize=onecolumn_figsize)
    axes = fig.add_subplot(111)

    if include_fast_tree and ft_corr is not None:
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

    if figname is not None:
        fig.savefig("{}.svg".format(figname))
        fig.savefig("{}.png".format(figname))
        fig.savefig("{}.pdf".format(figname))

if __name__ == '__main__':

    """
    Only trees generated using the following  simulation parameters will be used:
    """
    Mu=['0.0001'] # mutation rates
    Ts=['50'] # Frequency of sampling  during the simulations (see generate_sim_data script for details)


    """
    Should the FastTree accuracy be included in the plot
    """
    INCLUDE_FAST_TREE=True

    """
    should save figure
    """
    SAVE_FIG=False

    ##
    ## Input directories
    ##
    beast_trees_dir = "./simulated_data/2017-06-28_beast"
    other_trees_dir = "./simulated_data/dataset_06-28"

    ##
    ## Scan directory, read all trees, correlate similar branches
    ##
    tt_corr, ft_corr, bt_corr = correlation_dataset(other_trees_dir, beast_trees_dir, Mu=Mu, Ts=Ts)

    ##
    ## Plot the results
    ##
    plot_correlation(tt_corr, ft_corr, bt_corr,
        figname="./figs/simdata_BL_Corr" if SAVE_FIG else None
        include_fast_tree=INCLUDE_FAST_TREE)





