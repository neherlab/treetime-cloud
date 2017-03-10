#!/usr/bin/env python

from beast import run_beast
import os,sys
import datetime
import subprocess
import re
if __name__ == '__main__':

    tree = sys.argv[1]
    out_dir = sys.argv[2]
    print ("beast_run.py called for tree: " + tree)

    run_beast(tree, out_dir)
