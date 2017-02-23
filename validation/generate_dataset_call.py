#!/usr/bin/env python
"""
Actual script to perform LSd/TreeTime computations for a given set of parameters
"""
import os, sys
import generate_dataset as gds

if  __name__ == '__main__':

    sys.stderr.write ("  ".join(sys.argv) + "\n")

    L = int(float(sys.argv[1]))
    N = int(sys.argv[2])
    SAMPLE_VOL = int (sys.argv[3])
    SAMPLE_NUM = int(sys.argv[4])
    SAMPLE_FREQ = int(sys.argv[5])
    MU =  float(sys.argv[6])
    res_dir = sys.argv[7]
    suffix = sys.argv[8]
    outfile_prefix = sys.argv[9]

    sys.stderr.write  ("Importing GenDataSet... \n")

    sys.stderr.write (" Starting TreeTime run...\n ")

    # run evolution simulation, produce basic tree and alignment
    basename = gds.run_ffpopsim_simulation(L, N, SAMPLE_VOL, SAMPLE_NUM, SAMPLE_FREQ, MU, res_dir, suffix, optimize_branch_len=True)
    gds.reconstruct_fasttree(basename)

    # run treetime for original tree:
    outfile = outfile_prefix + "_treetime_res.csv"
    gds.run_treetime(basename, outfile, fasttree=False, failed=None)

    # run treetime for original tree:
    outfile = outfile_prefix + "_treetime_fasttree_res.csv"
    gds.run_treetime(basename, outfile, fasttree=True, failed=None)

    # run LSD for original tree:
    #directory to store all other LSD results
    if not os.path.exists(outfile_prefix + "_lsd"):
        os.mkdir(outfile_prefix + "_lsd")

    # file to store formatted results of LSD run
    outfile = outfile_prefix + "_lsd_res.csv"
    treefile = basename + ".opt.nwk"
    lsd_res_file = os.path.join(outfile_prefix + "_lsd", os.path.split(basename)[-1])
    gds.run_lsd(treefile, basename+".lsd_dates.txt", lsd_res_file, outfile)

    # run LSD for fast-tree tree:


    outfile = outfile_prefix + "_lsd_fasttree_res.csv"
    treefile = basename + ".ft.nwk"
    lsd_res_file = os.path.join(outfile_prefix + "_lsd", os.path.split(basename)[-1] + "_fasttree")
    gds.run_lsd(treefile, basename+"lsd_dates.txt", lsd_res_file, outfile)

    print ("Done")

