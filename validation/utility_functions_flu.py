#!/usr/bin/env python
"""
This module defines functions to facilitate operations with data specific
to Flu trees and alignments.
"""

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
import StringIO
import treetime
from utility_functions_general import remove_polytomies

def date_from_seq_name(name):
    """
    Parse flu sequence name to the date in numeric format (YYYY.F)

    Args:

     - name(str): name of the flu sequence.

    Returns:

     -  sequence sampling date if succeeded to parse. None otherwise.
    """
    def str2date_time(instr):
        """
        Convert input string to datetime object.
        Args:
         - instr (str): input string. Accepts one of the formats:
         {MM.DD.YYYY, MM.YYYY, MM/DD/YYYY, MM/YYYY, YYYY}.

        Returns:
         - date (datetime.datetime): parsed date object. If the parsing failed,
         None is returned
        """

        instr = instr.replace('/', '.')
        # import ipdb; ipdb.set_trace()
        try:
            date = datetime.datetime.strptime(instr, "%m.%d.%Y")
        except ValueError:
            date = None
        if date is not None:
            return date

        try:
            date = datetime.datetime.strptime(instr, "%m.%Y")
        except ValueError:
            date = None

        if date is not None:
            return date

        try:
            date = datetime.datetime.strptime(instr, "%Y")
        except ValueError:
            date = None

        return date

    try:
        date = str2date_time(name.split('|')[3].strip())
        return date.year + (date - datetime.datetime(date.year, 1, 1)).days / 365.25
    except:
        return None

def dates_from_flu_tree(tree):
    """
    Iterate over the Flu tree, parse each leaf name and return dates for the
    leaves as dictionary.

    Args:

     - tree(str or Biopython tree): Flu tree

    Returns:

     - dates(dict): dictionary of dates in format {seq_name: numdate}. Only the
     entries which were parsed successfully are included.
    """

    if isinstance(tree, str):
        tree = Phylo.read(tree, 'newick')

    dates = {k.name:date_from_seq_name(k.name) for k in tree.get_terminals()
                if date_from_seq_name(k.name) is not None}
    return dates

def subtree_with_same_root(tree, Nleaves, outfile, optimize=True):
    """
    Sample subtree of the given tree so that the root of the subtree is that of
    the original tree.

    Args:

     - tree(str or Biopython tree): initial tree

     - Nleaves(int): number of leaves in the target subtree

     - outfile(str): path to save the resulting subtree

     optimize(bool): perform branch length optimization for the subtree?

    Returns:
     - tree(Biopython tree): the subtree
    """
    if isinstance(tree, str):
        treecopy = Phylo.read(tree, 'newick')
    else:
        treecopy = copy.deepcopy(tree)

    remove_polytomies(treecopy)
    assert(len(treecopy.root.clades) == 2)

    tot_terminals = treecopy.count_terminals()

    # sample to the left of the root
    left = treecopy.root.clades[0]
    n_left = left.count_terminals()
    right = treecopy.root.clades[1]
    n_right = right.count_terminals()


    n_left_sampled = np.min((n_left, Nleaves * n_left / (n_left + n_right)))
    n_left_sampled = np.max((n_left_sampled, 5))  # make sure we have at least one
    left_terminals = left.get_terminals()
    left_sample_idx = np.random.choice(np.arange(len(left_terminals)), size=n_left_sampled, replace=False)
    left_sample = [left_terminals[i] for i in left_sample_idx]

    # sample to the right of the root
    n_right_sampled = np.min((n_right, Nleaves * n_right / (n_left + n_right)))
    n_right_sampled = np.max((n_right_sampled, 5))  # make sure we have at least one
    right_terminals = right.get_terminals()
    right_sample_idx = np.random.choice(np.arange(len(right_terminals)), size=n_right_sampled, replace=False)
    right_sample = [right_terminals[i] for i in right_sample_idx]

    for leaf in treecopy.get_terminals():
        if leaf not in right_sample and leaf not in left_sample:
            treecopy.prune(leaf)
        else:
            pass
            #print ("leaving leaf {} in the tree".format(leaf.name))

    if optimize:
        import treetime
        dates = dates_from_flu_tree(treecopy)
        aln = './resources/flu_H3N2/H3N2_HA_2011_2013.fasta'
        tt = treetime.TreeAnc(tree=treecopy, aln=aln,gtr='Jukes-Cantor')
        tt.optimize_seq_and_branch_len(prune_short=False)
        Phylo.write(tt.tree, outfile, 'newick')
        return tt.tree
    else:
        Phylo.write(treecopy, outfile, 'newick')
        return treecopy

def subtree_year_vol(tree, N_per_year, outfile):
    """
    Sample subtree of the given tree with equal number of samples per year.

    Note:

     - if there are not enough leaves sampled at a given year, all leaves for this
     year will be included in the subtree.

    Args:

     - tree(str or Biopython object): Initial tree

     - N_per_year(int): number of samples per year.

     - outfile (str): path to save the subtree

    Returns:
     - tree(Biopython tree): the subtree
    """

    if isinstance(tree, str):
        treecopy = Phylo.read(tree, 'newick')
    else:
        treecopy = copy.deepcopy(tree)

    remove_polytomies(treecopy)

    dates = dates_from_flu_tree(treecopy)
    sample = []

    cntr = Counter(map (int, dates.values()))
    years = cntr.keys()
    min_year = np.min(years)
    for year in years:
        all_names = [k for k in dates if int(dates[k]) == year]
        if len(all_names) <= N_per_year or year == min_year:
            sample += all_names
        else:
            sample += list(np.random.choice(all_names, size=N_per_year, replace=False))


    for leaf in treecopy.get_terminals():
        if leaf.name not in sample:
            treecopy.prune(leaf)
        else:
            pass
            #print ("leaving leaf {} in the tree".format(leaf.name))

    Phylo.write(treecopy, outfile, 'newick')
    return treecopy

def create_LSD_dates_file_from_flu_tree(tree, outfile):
    """
    Parse dates from the flu tree and write to the file in the LSD format.

    Args:

     - tree(str or Biopython object): Initial tree

     - outfile(str): path to save the LSD dates file.

    Returns:

     - dates(dict): dates parsed from the tree as dictionary.
    """

    dates = dates_from_flu_tree(tree)

    with open(outfile, 'w') as df:
        df.write(str(len(dates)) + "\n")
        df.write("\n".join([str(k) + "\t" + str(dates[k]) for k in dates]))
    return dates

def create_treetime_with_missing_dates(alnfile, treefile, dates_known_fraction=1.0):
    """
    Create TreeTime object with fraction of leaves having no sampling dates.
    The leaves to earse sampling dates are chosen randomly.

    Args:

     - alnfile(str): path to the flu alignment

     - treefiule(str): path to the Flu newixk tree

     - dates_known_fraction(float): fraction of leaves, which should have
     sampling date information.

    """
    aln = AlignIO.read(alnfile, 'fasta')
    tt = Phylo.read(treefile, 'newick')

    dates = {k.name: date_from_seq_name(k.name) for k in aln}

    # randomly choose the dates so that only the  known_ratio number of dates is known
    if dates_known_fraction != 1.0:
        assert(dates_known_fraction > 0 and dates_known_fraction < 1.0)
        knonw_keys = np.random.choice(dates.keys(), size=int (len(dates) * dates_known_fraction), replace=False)
        dates = {k : dates[k] for k in knonw_keys}

    myTree = treetime.TreeTime(gtr='Jukes-Cantor', tree = tt,
            aln = aln, verbose = 4, dates = dates, debug=False)

    myTree.optimize_seq_and_branch_len(reuse_branch_len=True, prune_short=True, max_iter=5, infer_gtr=False)

    return myTree

def create_subtree(tree, n_seqs, out_file, st_type='equal_sampling'):
    """
    Args:
     - tree(filename or Biopython tree): original tree

     - n_seqs: number of leaves in the resulting subtree

     - out_file: output locaton to store the resulting subtree

     - st_type: type of the subtree generation algorithm. Available types:
         - random: just choose n_leaves randomly
         - equal_sampling: choose equal leaves each year (if possible)
         - preserve_root: sample from right and left subtrees of the tree root.
         The root of the resulting subtree is therefore the same as of the original tree

    """
    if isinstance(tree, str):
        tree = Phylo.read(tree, 'newick')

    pass

if __name__ == '__main__':
    pass




