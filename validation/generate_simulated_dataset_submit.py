#!/usr/bin/env python
"""
Small script file to generate full dataset for the LSD/TreeTime validation on the
generated data. In order to perform parallel computations on a cluster, the script is
split into two parts. This part defines the input parameters, and for each set of the
parameters, calls the 'generate_dataset_call' script as a separate process.

Edit the subprocess call part according to the real configuration of the
computational facilities used. As an example, there is a direct subprocess call
commented, which can be used to perform computations on a local machine.
"""
import sys, os
import subprocess as sp
sys.path.append("./")

if __name__ == '__main__':

    CLUSTER = True

    # Directory to store results
    res_dir = "./simulated_data/dataset"
    # File prefix to store the formatted output for further processing/comparison
    outfile = "./simulated_data/2017-05-16"

    # FFPopSim simulation parameters
    L = 1e4
    N = 100
    SAMPLE_VOL = 10
    SAMPLE_NUM = 20
    SAMPLE_FREQS = [10, 20, 50] # T/N = [2, 4, 10]
    MUS = [7e-6, 1e-5, 2e-5, 5e-5, 7e-5, 1e-4, 2e-4, 5e-4, 7e-4, 1e-3, 2e-3]
    # MUS= [2e-4, 5e-4, 7e-4, 1e-3, 2e-3]
    N_POINTS = 20
    N_0 = 0

    # run treetime in-place:
    Ncalls = 0
    for MU in MUS:
        for SAMPLE_FREQ in SAMPLE_FREQS:
            for i in xrange(N_0, N_0 + N_POINTS):
                suffix = str(i)

                Ncalls += 1
                if Ncalls > 3000:
                    print ("Number of jobs exceeded")
                    break

                if CLUSTER:
                    call = ['qsub', '-cwd', '-b','y',
                           '-l', 'h_rt=23:59:0', # BEAST might run long
                            #'-o', './stdout.txt',
                              #'-e', './stderr.txt',
                            '-l', 'h_vmem=50G', # BEAST requires A LOT
                             './generate_simulated_dataset_run.py']
                else:
                    call = ['./generate_simulated_dataset_run.py']

                # run computations on a cluster

                arguments = [str(L),
                            str(N),
                            str(SAMPLE_VOL),
                            str(SAMPLE_NUM),
                            str(SAMPLE_FREQ),
                            str(MU),
                            res_dir,
                            suffix,
                            outfile]
                call.extend(arguments)

                sp.call(call)
