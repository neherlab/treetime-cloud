#!/usr/bin/env python
import treetime
import numpy as np
import os,sys
import datetime
import subprocess
import re

import utility_functions_flu as flu_utils
import utility_functions_general as gen_utils

aln_name = "./resources/flu_H3N2/H3N2_HA_2011_2013.fasta"
tree_name = "./resources/flu_H3N2/H3N2_HA_2011_2013.nwk"

run_treetime = True
run_LSD = True
run_Beast = True

if __name__ == "__main__":

    N_leaves = int(sys.argv[1])
    subtrees_dir = sys.argv[2]
    subtree_fname_suffix = sys.argv[3]
    treetime_res_file = sys.argv[4]

    lsd_res_file = sys.argv[5]
    lsd_outdir = sys.argv[6]

    if len(sys.argv) > 7:
        lsd_params = sys.argv[7].split("|")
    else:
        lsd_params = ['-c', '-r', 'a', '-v']


    if not os.path.exists(subtrees_dir):
        os.mkdir(subtrees_dir)

    #  Sample subtree
    subtree_fname_format = "H3N2_HA_2011_2013_{}_{}.nwk".format(N_leaves, subtree_fname_suffix)
    subtree_filename = os.path.join(subtrees_dir, subtree_fname_format)
    tree = flu_utils.subtree_with_same_root(tree_name, N_leaves, subtree_filename)
    N_leaves = tree.count_terminals()

    #  CALL TREETIME
    if run_treetime:
        dates = flu_utils.dates_from_flu_tree(tree)
        myTree = treetime.TreeTime(gtr='Jukes-Cantor',
            tree=tree, aln=aln_name, dates=dates,
            debug=False, verbose=4)
        myTree.optimize_seq_and_branch_len(reuse_branch_len=True, prune_short=True, max_iter=5, infer_gtr=False)
        start = datetime.datetime.now()
        myTree.run(root='best', relaxed_clock=False, max_iter=3, resolve_polytomies=True, do_marginal=False)
        end = datetime.datetime.now()
        with open(treetime_res_file, 'a') as of:
            of.write("{},{},{},{},{},{},{}\n".format(
                subtree_filename,
                str(N_leaves),
                str(myTree.tree.root.numdate),
                str(myTree.date2dist.slope),
                str(myTree.date2dist.r_val),
                str(gen_utils.internal_regress(myTree)),
                str((end-start).total_seconds())    ))
        print ("TreeTime done!")
    else:
        print ("Skip TreeTime run")


    if run_LSD:
        #  run LSD for the subtree:
        if not os.path.exists(lsd_outdir):
            os.mkdir(lsd_outdir)
        lsd_outfile = os.path.join(lsd_outdir, subtree_fname_format.replace(".nwk", ".txt"))
        datesfile = os.path.join(lsd_outdir, subtree_fname_format.replace(".nwk", ".lsd_dates.txt"))
        flu_utils.create_LSD_dates_file_from_flu_tree(subtree_filename, datesfile)
        runtime = gen_utils.run_LSD(subtree_filename, datesfile, lsd_outfile, lsd_params)
        #  parse LSD results
        tmrca, mu, objective = gen_utils.parse_lsd_output(lsd_outfile)
        if float(mu) <= 0:
            return 
            
        with open(lsd_res_file, "a") as of:
            of.write(",".join([subtree_filename, str(N_leaves), tmrca, mu, runtime, objective]))
            of.write("\n")

        print ("LSD Done!")
    else:
        print ("Skip LSD run")

    # run BEAST


