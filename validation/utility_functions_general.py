import numpy as np
from Bio import AlignIO, Phylo
from Bio.Align import  MultipleSeqAlignment
import random
import subprocess
import datetime
import os, copy
#import matplotlib.pyplot as plt
from scipy.stats import linregress
from collections import Counter
import StringIO

#plt.ion()
#plt.show()

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

def parse_lsd_output(lsd_outfile):

    # parse LSD results
    try:
        with open (lsd_outfile, 'r') as inf:
            ss = inf.readlines()
    except:

        folder, filename = os.path.split(lsd_outfile)
        files = [k for k in os.listdir(folder) if k.startswith(filename)]
        if len(files) == 0:
            return -1, -1, 0.0
        else:
            file = os.path.join(folder, files[0])
            with open(file, 'r') as inf:
                ss = inf.readlines()

    mystrs = [i for i in ss if 'tMRCA' in i]
    if len(mystrs) == 0:
        return -1, -1, 0.0

    mystr = mystrs[0]
    tmrca = (mystr.split(" ")[5])[:-1]
    mu = (mystr.split(" ")[3])[:-1]

    if len (mystr.split(" ")) == 8:
        objective = (mystr.split(" ")[7])[:-1]
    else:
        objective =  "0.0" # (mystr.split(" ")[7])

    return tmrca, mu, objective

def run_LSD(tree_filename, dates_filename, outfile, lsd_params = ['-r','a','-c','-v']):

    from external_binaries import LSD_BIN
    call = [LSD_BIN, '-i', tree_filename, '-d', dates_filename, '-o', outfile]
    call.extend(lsd_params)

    start = datetime.datetime.now()
    subprocess.call(call)
    end = datetime.datetime.now()
    runtime = str((end-start).total_seconds())
    return runtime

if __name__ == '__main__':

    pass
