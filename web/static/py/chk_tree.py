#!/home/pavel/.conda2/bin/python
import sys
from Bio import Phylo
import json 


def write_ok(out, t):
    data = {}
    data['status'] = 'OK'
    data['N_leaves'] = len(t.get_terminals())
    with open(out, 'w') as outf:
        json.dump(data, outf)

def write_err(out, string):
    data = {}
    data['status'] = 'error'
    data['desc'] = string
    with open(out, 'w') as outf:
        json.dump(data, outf)


if __name__ == '__main__':

    dir_ = sys.argv[1]
    out = sys.argv[2]

    print ("Python tree checking file: " + dir_)
    
    try:
        t = Phylo.read(dir_, 'newick')
        
        if len(t.get_terminals()) == 0 : 
            write_err(out, "Tree is empty!")
        else:
            print ("reading tree success..." + str(t.get_terminals()))
            sys.stdout.flush()
            write_ok(out, t)

    except:
        write_err(out, "Cannot read tree file: bad format")
        print("{'result':'error.'}")
        sys.stdout.flush()
