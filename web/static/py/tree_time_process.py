from __future__ import print_function, division
from treetime_web import TreeTimeWeb
import os, sys
dirname = (os.path.dirname(__file__))

def build_tree(root):

    aln_filename = os.path.join(root, "in_aln.fasta")
    tree_filename = os.path.join(root, "in_tree.nwk")
    fast_tree_bin = os.path.join(dirname, "fasttree")
    if not os.path.exists(fast_tree_bin):
        raise (RuntimeError("Cannot find FastTree binary."))

    call = [fast_tree_bin, '-nt','-quiet', aln_filename, ' > ', tree_filename]

    res = os.system(' '.join(call))
    if res != 0:
        raise RuntimeError("FastTree: Exception caught while building the NJ tree")

    return

def process(root):

    myTree = TreeTimeWeb(root, build_tree)
    myTree.run()


if __name__=="__main__":

    root = "./"
    process(root)



