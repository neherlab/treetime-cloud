import numpy as np
from Bio import AlignIO, Phylo
from Bio.Align import  MultipleSeqAlignment
import random
import subprocess
import datetime
import os, copy
import matplotlib.pyplot as plt
from scipy.stats import linregress
plt.ion()
plt.show()

def date_from_seq_name(name):
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

    date = str2date_time(name.split('|')[2].strip())

    return date.year + (date - datetime.datetime(date.year, 1, 1)).days / 365.25

def remove_polytomies(tree):
    """
    Scan tree, and remove the polytomies (if any) by random merging
    Returns: tree without polytomies
    """
    for clade in tree.find_clades():
        #if not hasattr(clade, "name") or clade.name is None:
        #    clade.name = "None"
        if len(clade.clades) < 3:
            continue
        clades = clade.clades
        while len(clades) > 2:
            c1 = clades.pop()
            c2 = clades.pop()
            new_node = Phylo.BaseTree.Clade()
            new_node.name="None"
            new_node.branch_length = 1e-5
            new_node.clades = [c1,c2]
            clades.append(new_node)
        clade.clades = clades

    return tree

def subtree_with_same_root(tree, Nleaves, outfile):

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
    left_sample_idx = np.random.choice(np.arange(len(left_terminals)), size=n_left_sampled)
    left_sample = [left_terminals[i] for i in left_sample_idx]

    # sample to the right of the root
    n_right_sampled = np.min((n_right, Nleaves * n_right / (n_left + n_right)))
    n_right_sampled = np.max((n_right_sampled, 5))  # make sure we have at least one
    right_terminals = right.get_terminals()
    right_sample_idx = np.random.choice(np.arange(len(right_terminals)), size=n_right_sampled)
    right_sample = [right_terminals[i] for i in right_sample_idx]

    for leaf in treecopy.get_terminals():
        if leaf not in right_sample and leaf not in left_sample:
            treecopy.prune(leaf)
        else:
            pass
            #print ("leaving leaf {} in the tree".format(leaf.name))

    Phylo.write(treecopy, outfile, 'newick')
    return treecopy

def dates_from_flu_tree(tree):

    if isinstance(tree, str):
        tree = Phylo.read(tree, 'newick')

    dates = {k.name:date_from_seq_name(k.name) for k in tree.get_terminals()
                if date_from_seq_name(k.name) is not None}
    return dates

def LSD_dates_file_from_tree(tree, outfile):

    dates = dates_from_flu_tree(tree)

    with open(outfile, 'w') as df:
        df.write(str(len(dates)) + "\n")
        df.write("\n".join([str(k) + "\t" + str(dates[k]) for k in dates]))
    return dates

def internal_regress(myTree):
    resarr = []
    for node in myTree.tree.get_nonterminals():
        try:
            resarr.append((node.numdate, node.dist2root))
        except:
            continue
    resarr = np.array(resarr)
    if resarr.shape[0] == 0:
        return 0.
    else:
        return linregress(resarr[:, 0], resarr[:, 1]).rvalue**2

if __name__ == '__main__':

    tree = subtree_with_same_root('./H3N2_HA_1980_2015_NA.nwk', 1000, './H3N2_HA_1980_2015_NA_100.nwk')




