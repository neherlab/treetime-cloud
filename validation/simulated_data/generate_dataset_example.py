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

    # Directory to store results
    res_dir = "./accuracy_5"
    # File prefix to store the formatted output for further processing/comparison
    outfile = "./accuracy_5_"

    # FFPopSim simulation parameters
    L = 1e4
    N = 100
    SAMPLE_VOL = 10
    SAMPLE_NUM = 20
    SAMPLE_FREQS = [10, 20, 50, 75, 100] #[1, 2, 3, 5, 10, 20, 50, 75, 100, 150, 200]
    MUS = [1e-5, 2e-5, 5e-5, 8e-5, 1e-4, 2e-4, 5e-4, 7e-4, 1e-3, 2e-3, 4e-3]
    N_POINTS = 20
    N_0 = 30

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

                # direct invokation of a subprocess
                #call = ['./generate_dataset_call.py']

                # run computations on a cluster
                call = ['qsub', '-cwd', '-b','y',
                       '-l', 'h_rt=0:59:0',
                        #'-o', './stdout.txt',
                        #'-e', './stderr.txt',
                        '-l', 'h_vmem=3G',
                         './generate_dataset_call.py']

                arguments = [str(L), str(N), str(SAMPLE_VOL), str(SAMPLE_NUM),
                             str(SAMPLE_FREQ), str(MU), res_dir, suffix, outfile]
                call.extend(arguments)

                sp.call(call)


