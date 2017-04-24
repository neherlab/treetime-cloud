#!/usr/bin/env python

import treetime
import numpy as np
import os,sys
import analysis
import datetime
import subprocess
import re

aln_name = "./H3N2_HA_1980_2015_NA.fasta"
LSD_BIN = "/ebio/ag-neher/share/users/psagulenko/programs/LSD/lsd-0.2/bin/lsd"

if __name__ == "__main__":


    work_dir = sys.argv[1]
    N_per_year = int(sys.argv[2])
    treename = sys.argv[3]
    res_file = sys.argv[4]
    lsd_res_file = sys.argv[5]
    filename_suffix = sys.argv[6]

    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
    fname_format = "H3N2_HA_1980_2015_NA_{}_{}.nwk".format(N_per_year, filename_suffix)
    filename = os.path.join(work_dir, fname_format)
    #  Sample subtree
    tree = analysis.subtree_year_vol(treename, N_per_year, filename)
    N_leaves = tree.count_terminals()

    failed = []
    try:

        ## CALL TREETIME
        dates = analysis.dates_from_flu_tree(tree)
        myTree = treetime.TreeTime(gtr='Jukes-Cantor',
            tree=tree, aln=aln_name, dates=dates,
            debug=False, verbose=4)

        myTree.optimize_seq_and_branch_len(self,reuse_branch_len=True, prune_short=True, max_iter=5, infer_gtr=False)

        start = datetime.datetime.now()
        #
        myTree.run(root='best', relaxed_clock=False, max_iter=1, resolve_polytomies=True, do_marginal=False)
        #
        end = datetime.datetime.now()
        with open(res_file, 'a') as of:
            of.write("{},{},{},{},{},{},{}\n".format(
                filename,
                str(N_leaves),
                str(myTree.tree.root.numdate),
                str(myTree.date2dist.slope),
                str(myTree.date2dist.r_val),
                str(analysis.internal_regress(myTree)),
                str((end-start).total_seconds()) ))

    except:
        if failed is not None:
            failed.append(basename)

    ##  CALL LSD
    lsd_outdir = "./LSD"
    # run LSD for original tree:
    #directory to store all other LSD results
    if not os.path.exists(lsd_outdir):
        os.mkdir(lsd_outdir)

    lsd_outfile = os.path.join(lsd_outdir, fname_format.replace(".nwk", ".txt"))
    lsd_tree = filename.replace(".nwk", ".opt.nwk")
    datesfile = os.path.join(lsd_outdir, fname_format.replace(".nwk", ".lsd_dates.txt"))
    analysis.LSD_dates_file_from_tree(filename, datesfile)
    #import ipdb; ipdb.set_trace()
    call = [LSD_BIN, '-i', filename, '-d', datesfile, '-o', lsd_outfile,
                '-c',
                '-r', 'a',
                '-v']

    start = datetime.datetime.now()
    subprocess.call(call)
    end = datetime.datetime.now()

    # parse LSD results
    try:
        with open (lsd_outfile, 'r') as inf:
            ss = inf.readlines()
    except:
        with open (lsd_outfile + 'q', 'r') as inf:
            ss = inf.readlines()


    mystr = [i for i in ss if 'tMRCA' in i][0]
    tmrca = (mystr.split(" ")[5])[:-1]
    mu = (mystr.split(" ")[3])[:-1]
    runtime=str((end-start).total_seconds())

    if len (mystr.split(" ")) == 8:
        objective = (mystr.split(" ")[7])[:-1]
    else:
        objective =  "0.0" # (mystr.split(" ")[7])

    with open(lsd_res_file, "a") as of:
        of.write(",".join([filename, str(N_leaves), tmrca, mu, runtime, objective]))
        of.write("\n")

