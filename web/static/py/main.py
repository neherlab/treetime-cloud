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
            'treetime':tt,
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
            'treetime':tt,
            'slope': slope,
            'branch_penalty':bp,
            'do_reroot':do_reroot            
        }))

    # 6. run treetime
    steps.append(create_step(
        "Run TreeTime",
        run_treetime,
        {
            'treetime':tt
        }))
    # 7. relaxd clock
    if is_true(s, 'doRelaxedClock'):
        steps.append(create_step(
            "Calc relaxed molecular clock",
            relaxed_clock,
            {
                'treetime':tt,
                alpha: get_value(s, 'doRelaxedClock', 'alpha'),
                beta: get_value(s, 'doRelaxedClock', 'beta'),
            }))

    return steps


# steps - unit functions
def basic_setup(state, root):
    state['status'] = 'Done'
    return True
    pass

def build_tree(state, aln):
    state['status'] = 'Done'
    print ('aln: ' + aln)
    return True
    pass

def init_treetime(state, root): #  read tree, aln, met, return treetime object
    state['status'] = 'Done'
    return True
    pass

def anc(state, treetime, reuse_branch=False):
    state['status'] = 'Done'
    return True
    pass

def temporal_constraints(state, treetime, slope=None, branch_penalty=0.0, do_reroot=False):
    state['status'] = 'Done'
    return True
    pass

def run_treetime(state, treetime):
    state['status'] = 'Done'
    return True
    pass

def resolve_poly(state, treetime):
    state['status'] = 'Done'
    return True
    pass

def relaxed_clock(state, treetime, alpha, beta):
    state['status'] = 'Done'
    return True
    pass

def coalescent_models(state, treetime, Tc):
    state['status'] = 'Done'
    return True
    pass

def root_variance(state, treetime):
    state['status'] = 'Done'
    return True
    pass

def do_gtr_calc(state, treetime):
    state['status'] = 'Done'
    return True
    pass

if __name__ == '__main__':
    

    print("Passed")
    sys.stdout.flush()

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
    for step in steps:
        settings = step['settings']
        call = step['callable']
        res = call(state[idx], **settings)
        write_json(ostate, state)
        if not res:
            break
        idx += 1

    sys.exit(1)

    #  files 
    session_state = root + "session_state.json"
    nwk = root + "in_tree.nwk"
    aln = root + "in_aln.fasta"
    meta = root + "in_meta.csv"
    
    with open(session_state, 'r') as inf:    
        ss_json = json.load(inf)

    #  steps 
    ss_json['steps'] = {
        '0':{'status': 'Todo', 'name': 'Basic Setup'},
        #'1':{'status': 'Todo', 'name': 'Build Tree'},
        '2':{'status': 'Todo', 'name': 'Initialize TreeTime objects'},
        '3':{'status': 'Todo', 'name': 'Ancestral reconstruction'},
        '4':{'status': 'Todo', 'name': 'Set temporal constraints'},
        '6':{'status': 'Todo', 'name': 'Run TreeTime'}
    }

    steps = ss_json['steps']

    should_build_tree = ss_json['run']['should_build_tree'] == 'true' or ss_json['run']['should_build_tree'] == True
    should_use_branch_penalty = ss_json['run']['should_use_branch_penalty'] == 'true' or ss_json['run']['should_use_branch_penalty'] == True
    do_root_calc = ss_json['run']['do_root_calc'] == 'true' or ss_json['run']['do_root_calc'] == True
    should_use_slope =  ss_json['run']['should_use_slope'] == True or ss_json['run']['should_use_slope'] == 'true' 
    # add optional steps to the workflow:
    if should_build_tree:
        steps['1'] = {'status': 'Todo', 'name': 'Build Tree'},
    
    if do_root_calc:
        steps['5'] = {'status': 'Todo', 'name' : 'Optimize root position'}

    

    # save json
    write_json(session_state, ss_json)

    # initialization for the run parameters
    try:
        sys.stdout.write ("reading run_params parameters...")
        run_params = ss_json['run']
        # todo check all files

    except:
        steps['0']['status'] = "Failed"
        write_json(session_state, ss_json)
        sys.exit(0)
    
    # reading run_params parameters
    steps['0']['status'] = "Done"

    # build tree
    try:
        if should_build_tree: 
            steps['1']['status'] = 'Progress'
            write_json(session_state, ss_json)
            # TODO build tree with fast tree
        
            steps['1']['status'] ='Done'
            write_json(session_state, ss_json)
            
    
    except:

        sys.stdout.write ('Basic Setup failed!')
        steps['1']['status'] = 'Failed'
        steps['1']['message'] = 'Internal server error'
        write_json(session_state, ss_json)
        sys.exit(0)

    # tree file OK, tree has been built
    #  initialize basi treetime
    try:

        steps['2']['status'] = 'Progress'
        write_json(session_state, ss_json)
        
        if run_params['gtr'] not in gtrs:
            steps['2']['status'] = 'Failed'
            steps['2']['message'] = 'GTR model of type' + run_params['gtr'] + 'is not defined'
            write_json(session_state, ss_json)
            sys.exit(0)

        # create gtr model
        gtr = create_GTR(run_params['gtr'])
        t = treetime.treetime_from_newick(gtr, nwk)
        treetime.set_seqs_to_leaves(t, AlignIO.read(aln, 'fasta'))
        treetime.read_metadata(t, meta)

        steps['2']['status'] = 'Done'
        write_json(session_state, ss_json)

    except:

        steps['2']['status'] = 'Failed'
        steps['2']['message'] = 'Internal server error in reading input files'
        write_json(session_state, ss_json)
        sys.stdout.write ('Initialization of the  TreeTime parameters failed!')
        sys.exit(0)

    # initialization went successful - run optimization of the branch lengths
    try:
        
        steps['3']['status'] = 'Progress'
        write_json(session_state, ss_json)
        
        should_use_branch_len = run_params['should_use_branch_len']
        t.optimize_seq_and_branch_len(should_use_branch_len, True)
        
        steps['3']['status'] = 'Done'
        write_json(session_state, ss_json)
    
    except:
        
        steps['3']['status'] = 'Failed'
        steps['3']['message'] = 'TreeTime script error'
        write_json(session_state, ss_json)
        sys.stdout.write ('Ancestral reconstruction failed!')
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

        t.ml_t()

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
    

    


    

