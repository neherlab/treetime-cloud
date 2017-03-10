#!/usr/bin/env python

from generate_dataset import run_beast
import sys
if __name__ == '__main__':

    tree = sys.argv[1]
    out_dir = sys.argv[2]
    print ("beast_run.py called for tree: " + tree)
    assert(tree.endswith ('.opt.nwk'))
    run_beast(tree, out_dir)
