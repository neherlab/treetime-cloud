#!/home/pavel/.conda2/bin/python

from __future__ import print_function, division
import treetime
import numpy as np
import datetime
import os,sys,copy,json
from Bio import Phylo, AlignIO
import matplotlib.pyplot as plt
plt.ion()
plt.show()

gtrs = ['jc'] # available gtr 
tt = None # basic treetime object

def create_GTR(gtr):
    if gtr != 'jc': 
        raise NotImplementedError("")
    else: 
        return treetime.GTR.standard()

def write_json(session_state, ss_json):
    with open(session_state, 'w') as of:
            json.dump(ss_json, of)


def create_step(name, call, settings):
    return {"name":name,
    "callable":call,
    "settings":settings}
    pass


def load_settings(root):
    
    def is_true(s, entry):
        
        if entry not in s: 
            return False
        return s[entry] == 'true' or  s[entry] == 'True' or s[entry] is True

    def get_value(s,name,val):
        
        res = None
        try:
            if name in s:
                _should = s[name]
                _bool = is_true(_should, 'bool')
            if _bool and val in _should:
                res = float(_should[val])
        except:
            print ("Error in {%0} parsing! {%2} is NONE").format((name, val))
            res = None
        return res


    with open(os.path.join(root, 'settings.json')) as ff:
        s = json.load(ff)
    
    steps = []
    
    # 1. basic setup
    steps.append(create_step("Basic setup", basic_setup, {'root':root}))
    
    # 2. decide should we build tree
    if is_true(s, 'doBuildTree'):
        steps.append(create_step("Build initial tree", 
                                 build_tree, 
                                 {'aln':os.path.join(root,"in_aln.fasta")}))       

    # 3. init treetime objects
    steps.append(create_step("Initialize tree-time objects", 
                             init_treetime,
                             {
                                'root':root
                             }))
    
    # 4. ancestral state reconstruction:
    steps.append(create_step(
        "Ancestral state reconstruction", 
        anc, 
        {
            'reuse_branch': is_true(s, 'shouldReuseBranchLen')
        }))

    # 5. init temporal constraints
    slope = get_value(s, 'shouldUseSlope', 'value')
    bp = get_value(s, 'shouldUseBranchLenPenalty', 'value')
    do_reroot = is_true(s, 'doReroot')
    if do_reroot:
        name = "Optimizing root position, initializing temporal constraints"
    else:
        name = "Initializing temporal constraints"

    steps.append(create_step(
        name,
        temporal_constraints,
        {
            
            'slope': slope,
            'branch_penalty':bp,
            'do_reroot':do_reroot            
        }))

    # 6. run treetime
    steps.append(create_step(
        "Run TreeTime",
        run_treetime,
        {
            
        }))
    # 7. relaxd clock
    if is_true(s, 'doRelaxedClock'):
        steps.append(create_step(
            "Calc relaxed molecular clock",
            relaxed_clock,
            {
                
                alpha: get_value(s, 'doRelaxedClock', 'alpha'),
                beta: get_value(s, 'doRelaxedClock', 'beta'),
            }))

    return steps


# steps - unit functions
def basic_setup(tt, state, root):
    state['status'] = 'Done'
    return None, True
    pass

def build_tree(tt, state, aln):
    state['status'] = 'Done'
    print ('aln: ' + aln)
    return None, True
    pass

def init_treetime(tt, state, root): #  read tree, aln, met, return treetime object

    nwk  = os.path.join(root, "in_tree.nwk")
    aln  = os.path.join(root, "in_aln.fasta")
    meta = os.path.join(root, "in_meta.csv")
   
    try:
        # create gtr model
        gtr = treetime.GTR.standard()
        #  set global object for further use 
        tt = treetime.treetime_from_newick(gtr, nwk)
        treetime.set_seqs_to_leaves(tt, AlignIO.read(aln, 'fasta'))
        treetime.read_metadata(tt, meta)
        state['status'] = 'Done'
        return (tt, True)
    
    except:
        state['status'] = 'Error'
        return (None, False)
    

def anc(tt, state, reuse_branch=False):
    
#    try:
        tt.optimize_seq_and_branch_len(reuse_branch, True)
        state['status'] = 'Done'
        return (tt, True)
#    except:
        state['status'] = 'Error'
        return (None, False)

def temporal_constraints(tt, state, slope=None, branch_penalty=0.0, do_reroot=False):
    try:
        if branch_penalty != 0.0:
            treetime.ttconf.BRANCH_LEN_PENALTY = branch_penalty

        if do_reroot: # TODO slope also should be taken into account!
            tt.reroot_to_best_root()
        else:
            if slope == 0.0 : slope = None
            tt.init_date_constraints(slope=slope)


        state['status'] = 'Done'
        return tt, True
    except:
        return None, False

def run_treetime(tt, state):
    try:
        tt.ml_t()
        tt.tree.ladderize()
        state['status'] = 'Done'
        return tt, True
    except:
        state['status'] = 'Error'
        return None, False

def resolve_poly(tt, state):
    try:
        tt.resolve_polytomies()
        tt.tree.ladderize()
        state['status'] = 'Done'
        return tt, True
    except:
        state['status'] = 'Error'
        return None, False

def relaxed_clock(tt, state, alpha, beta):
    state['status'] = 'Done'
    return tt, True

def coalescent_models(tt, state, Tc):
    state['status'] = 'Done'
    return tt, True

def root_variance(tt, state):
    state['status'] = 'Done'
    return tt, True

def do_gtr_calc(tt, state):
    state['status'] = 'Done'
    return tt, True

if __name__ == '__main__':

    root = sys.argv[1]

    steps = load_settings(root)
    ostate = os.path.join(root, 'state.json')

    # write initial state 
    state = []
    for step in steps:
        name = step['name']
        state.append({
            'name': name,
            'status':'Todo', 
            'urgency': "", # warn or error
            'message':"" # in case of warn, error - message is important
        })
    write_json(ostate, state)


    # run the workflow
    idx = 0
    tt=None
    for step in steps:
        settings = step['settings']
        call = step['callable']
        tt, res = call(tt, state[idx], **settings)
        write_json(ostate, state)
        if not res:
            break
        idx += 1

    #  save files
    treetime.treetime_to_json(tt, "out_tree.json")
    treetime.tips_data_to_json(tt, "out_tips.json")
    treetime.root_lh_to_json(tt, "out_root_lh.json")

    sys.exit(0)


    #  Init date params, 
    #  if needed, do the following:
    #  correct root node
    #  use the slope value given 
    #  set branch length penalty 
    try:

        steps['4']['status']  = 'Progress'
        write_json(session_state, ss_json)

        if should_use_branch_penalty:
            treetime.treetime_conf.BRANCH_LEN_PENALTY = float(run_params['branch_penalty'])
        
        if not do_root_calc:
            if should_use_slope:
                slope = float(run_params['slope_value'])
            else:
                slope = None
            t.init_date_constraints(slope=slope)
        
        else: #  do root calc
            
            steps['5']['status'] =  'Progress'
            write_json(session_state, ss_json)
            
            if should_use_slope:
                slope = float(run_params['slope_value'])
            else:
                slope = None
            
            # infer better root
            t.reroot_to_best_root()
            #  TODO slope constraints

            steps['5']['status'] =  'Done'
            write_json(session_state, ss_json)
            

        steps['4']['status'] = 'Done'
        write_json(session_state, ss_json)
    
    except:
        
        steps['4']['status'] = 'Failed'
        steps['4']['message'] = 'TreeTime script error. Parameters are inconsistent.'
        write_json(session_state, ss_json)
        sys.stdout.write ('Init date params failed!')
        sys.exit(0)

    #  MAIN method
    try:
        steps['6']['status'] = 'Progress'
        write_json(session_state, ss_json)

        

        steps['6']['status'] =  'Done'
        write_json(session_state, ss_json)

    except:

        steps['6']['status'] = 'Failed'
        steps['6']['message'] =  'TreeTime script error. Parameters are inconsistent.'
        write_json(session_state, ss_json)
        sys.stdout.write ('TreeTime main failed!')
        sys.exit(0)

    #  TODO autocorr molecular clock
    #  TODO coalescence model
    #  TODO resolve poly


    print ("TreeTime python script Done!")
    

    


    

