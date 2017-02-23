#!/usr/bin/env python

import subprocess as sp
import numpy as np

if __name__ =="__main__":

    work_dir = "./subtrees"

    # leaves per year
    N_leaves_array = [1, 2, 3, 5, 10, 20, 30, 50, 100, 200] #[20, 50, 100, 150, 200, 250, 500, 1000, 1500, 2000, 3000, 4000, 5000]
    treename = "./H3N2_HA_1980_2015_NA.nwk"
    res_file = "H3N2_HA_1980_2015_RES_treetime.csv"
    lsd_res_file = "H3N2_HA_1980_2015_RES_lsd.csv"

    n_iter = 50
    for N_leaves in N_leaves_array:
        for iteration in np.arange(n_iter):

            #  uncomment to run benchmarks directly and synchronously
            #call = ['./benchmark_run.py']

            #  uncomment to run benchmarks in parallel on a cluster
            call = ['qsub', '-cwd', '-b','y',
                    '-l', 'h_rt=1:59:0',
                    #'-o', './stdout.txt',
                    #'-e', './stderr.txt',
                    '-l', 'h_vmem=3G',
                    './benchmark_run.py']

            arguments = [work_dir, str(N_leaves), treename, res_file, lsd_res_file, str(iteration)]
            call.extend(arguments)
            sp.call(call)