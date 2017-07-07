#!/usr/bin/env python
"""
Actual script to perform LSd/TreeTime computations for a given set of parameters
"""
import os, sys
import utility_functions_simulated_data as utils_sim

if  __name__ == '__main__':

    GENERATE_SIMULATED_DATA = False
    RUN_TREETIME = True
    RUN_LSD = True
    RUN_BEAST = True

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

    # run evolution simulation, produce basic tree and alignment
    if GENERATE_SIMULATED_DATA:
        print ("Running FFpopSim")
        basename = utils_sim.run_ffpopsim_simulation(L, N, SAMPLE_VOL, SAMPLE_NUM, SAMPLE_FREQ, MU, res_dir, suffix, optimize_branch_len=True)

        print ("Running FastTree reconstruction")
        utils_sim.reconstruct_fasttree(basename, optimize_branch_len=False)

    else:  # use previously generated data
        basename = "FFpopSim_L{}_N{}_Ns{}_Ts{}_Nv{}_Mu{}".format(str(L), str(N),
                 str(SAMPLE_NUM), str(SAMPLE_FREQ), str(SAMPLE_VOL), str(MU))

        basename = os.path.join(res_dir, basename)
        if suffix != "":
             basename = basename + "_" + suffix

    if RUN_TREETIME:
        # run treetime for original tree:
        print ("Running TreeTime, original tree")
        outfile = outfile_prefix + "_treetime_res.csv"
        utils_sim.run_treetime(basename, outfile, fasttree=False, failed=None, max_iter=3, use_input_branch_length=True)

        # run treetime for FastTree-reconstructed tree:
        print ("Running TreeTime, FastTree tree")
        outfile = outfile_prefix + "_treetime_fasttree_res.csv"
        utils_sim.run_treetime(basename, outfile, fasttree=True, failed=None, max_iter=3, use_input_branch_length=True)

        #utils_sim.run_treetime(basename, outfile_prefix + "_treetime_fasttree_res_use_input_branch_false.csv", fasttree=True, failed=None, max_iter=3, use_input_branch_length=False)

    if RUN_LSD:
        # run LSD for original tree:
        # directory to store all other LSD results
        try:
            if not os.path.exists(outfile_prefix + "_lsd"):
                os.mkdir(outfile_prefix + "_lsd")
        except:
            pass

        # file to store formatted results of LSD run
        print ("Running LSd, original tree")
        outfile = outfile_prefix + "_lsd_res.csv"
        treefile = basename + ".opt.nwk"
        lsd_res_file = os.path.join(outfile_prefix + "_lsd", os.path.split(basename)[-1])
        lsd_dates_file = os.path.join(outfile_prefix + "_lsd", os.path.split(basename)[-1] + ".lsd_dates.txt")
        utils_sim.run_lsd(treefile, lsd_dates_file, lsd_res_file, outfile)

        # run LSD for fast-tree tree:
        print ("Running LSD, FastTree tree")
        outfile = outfile_prefix + "_lsd_fasttree_res.csv"
        treefile = basename + ".ft.nwk"
        lsd_res_file = os.path.join(outfile_prefix + "_lsd", os.path.split(basename)[-1] + "_fasttree")
        utils_sim.run_lsd(treefile, lsd_dates_file, lsd_res_file, outfile)

    if RUN_BEAST:
        # run BEAST for the original tree:
        print ("Running BEAST, FastTree tree")
        outfile = outfile_prefix + "_beast_res.csv"
        utils_sim.run_beast(basename, out_dir=outfile_prefix+"_beast", res_file=outfile, fast_tree=True)
