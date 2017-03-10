#!/usr/bin/env python

import subprocess as sp
import os

if __name__ == '__main__':


    subtrees = './accuracy_5'
    out_dir = "./beast_out_cluster"


    trees = [os.path.join(subtrees, k) for k in os.listdir(subtrees) if k.endswith('_31.opt.nwk')]
    #trees = [os.path.join(subtrees, k) for k in os.listdir(subtrees) if k.endswith ('_1.nwk')]

    print ("\n".join(trees))

    for tree in trees:

        #  uncomment to run benchmarks directly and synchronously
        #call = ['./beast_run.py']
        #  uncomment to run benchmarks in parallel on a cluster
        call = ['qsub', '-cwd', '-b','y',
                '-l', 'h_rt=23:59:0',
                #'-o', './stdout.txt',
                #'-e', './stderr.txt',
                '-l', 'h_vmem=50G',
                './beast_run.py']


        arguments = [tree, out_dir]
        call.extend(arguments)
        sp.call(call)





