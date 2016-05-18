from __future__ import print_function, division
#import treetime # will be imported by the pipeline function
import numpy as np
import datetime
import os,sys,copy,json
from Bio import Phylo, AlignIO, Seq, SeqRecord, Align
from Bio.Alphabet import generic_dna
import matplotlib.pyplot as plt
import zipfile 
import StringIO

def _write_json(state_dic, outfile):
    with open(outfile, 'w') as of:
            json.dump(state_dic, of)

def create_initial_state(root, config_dic):
    """
    Create the list of steps to perform from the configuration supplied by the user
    """

    # create steps from settings
    state = {
        "state":"Running",
        "warn":[],
        "todo":[],
        "progress":"",
        "done":[],
    }
    if config_dic["do_build_tree"]==True:
        state["todo"].append("Build Phylogenetic tree")

    state["todo"].append("Basic TreeTime objects setup")

    if config_dic["reuse_branch_len"] == False:
        state["todo"].append("Reconstructing ancestral sequences, optimizing tree branches")

    if config_dic["reroot"] == True:
        # This function does the following:
        # 1. re-root 
        # 2. infer GTR model if needed
        # 3. ancestral inference
        # 4. init time constraints
        state["todo"].append("Search for the best root, infer GTR model, reconstruct ancestral sequences, initialize temporal constraints")        

    elif config_dic["infer_gtr"]==True: # 
        state["todo"].append("Infer GTR model, reconstruct ancestral sequences, initialize temporal constraints")
    else:
        state["todo"].append("Reconstruct ancestral sequences, initialize temporal constraints")

    state["todo"].append("Run TreeTime ML optimization")

    
    # some post-processing
    if config_dic["coalescent"]==True:
        state["todo"].append("Run coalescent model, resolve multiple mergers")
    elif config_dic["resolve_poly"]==True:
        state["todo"].append("Resolve multiple mergers")
    else:
        pass  # do nothing

    if config_dic["relax_mu"]==True:  # Relaxed molecular clock 
        state["todo"].append("Relax molecular clock")
    
    state["todo"].append("Save results")
    return state

def _update_session_state(root, state, state_fname="session_state.json", err=None, warn=None):
    
    if warn is not None:
        state["warn"].append(warn)
        return 

    if err is not None:
        state["error"] = err
    
    else: # normal run
        # update state dictionary
        if state["progress"] != "":
            state["done"].append(state["progress"])

        state["progress"] = ""
        
        if len(state["todo"]) < 1:
            state["state"]="complete"
        else:
            step_name = state["todo"][0]
            state["progress"] = step_name
            state["todo"].remove(step_name)

    _write_json(state, os.path.join(root, state_fname))
    print (state)
    # TODO save state to json file 

def build_tree(root):
    try:

        aln_filename = os.path.join(root, "in_aln.nwk")
        tree_filename = os.path.join(root, "in_tree.nwk")
    
        call = ['/ebio/ag-neher/share/programs/bin/fasttree', '-nt','-quiet',\
            aln_filename, ' > ', tree_filename]
            
        os.system(' '.join(call))
        return True
    
    except:

        return False

def process(root, config_dic, state_fname="session_state.json"):

    import treetime
    
    # some global configurations:
    if config_dic["use_mu"]==True:
        try:
            mu = float(config_dic["mu"])
        except:
            print ("Cannot set mutation rate from the string: %s. Data not understood" % config_dic["mu"])
            mu = None
    else:
        mu = None

    # specify filenames
    nwk = os.path.join(root, "in_tree.nwk")
    aln = os.path.join(root, "in_aln.fasta")
    meta = os.path.join(root, "in_meta.csv")

    #initialize the session state dictionary from the config_dic:
    state = create_initial_state(root, config_dic)
    
    err=None    
    if config_dic["do_build_tree"] == True:
        _update_session_state(root, state, state_fname=state_fname, err=err) 
        try:
            build_tree(root)
        except:
            err = "Error in TreeTime run: cannot build phylogenetic tree."
            _update_session_state(root, state, state_fname=state_fname, err=err)
            return

    _update_session_state(root, state, state_fname=state_fname, err=err)
    try:

        gtr = treetime.GTR.standard() # always start from standard J-C model
        tt = treetime.treetime_from_newick(gtr, nwk)  # load tree
        treetime.set_seqs_to_leaves(tt, AlignIO.read(aln, 'fasta'))  # assign sequences
        treetime.read_metadata(tt, meta)  # load meta data 
    
    except:
        err = "Cannot initialize TreeTime objects."
        _update_session_state(root, state, state_fname=state_fname, err=err)
        return


    ##  Pipeline
    if config_dic["reuse_branch_len"] == False:
        _update_session_state(root, state, state_fname=state_fname,  err=err)
        try:
            tt.optimize_seq_and_branch_len(False, True)
        except:
            err = "Branch lengths optimization failed."
            _update_session_state(root, state, state_fname=state_fname, err=err)
            return

    
    _update_session_state(root, state, state_fname=state_fname, err=err)
    try:
        if config_dic["reroot"] == True:
            # This function does the following:
            # 1. re-root 
            # 2. infer GTR model if needed
            # 3. ancestral inference
            # 4. init time constraints
            infer_gtr=config_dic["infer_gtr"]
            tt.reroot_to_best_root(infer_gtr=infer_gtr) # FIXME include mu slope here
        elif config_dic["infer_gtr"]==True: # 
            tt.infer_gtr()  # will make first shot in ancestral state reconstruction
            tt.init_date_constraints(slope=mu, ancestral_inference=True) 
        else:
            # just initialize  the objects needed to run the TreeTime-ML, 
            # and infer ancestral sequences
            tt.init_date_constraints(slope=mu, ancestral_inference=True)

    except:
        err = "Preparation of the TreeTime object failed."
        _update_session_state(root, state, state_fname=state_fname, err=err)
        return   

    # run TreeTime-ML optimization
    _update_session_state(root, state, state_fname=state_fname, err=err)
    try:
        tt.ml_t()
    except:
        err = "Error in the main TreeTime script. Cannot continue."
        _update_session_state(root, state, state_fname=state_fname, err=err)
        return   

    # some post-processing
    _update_session_state(root, state, state_fname=state_fname, err=err)
    try:
        if config_dic["coalescent"]==True:
            tt.coalescent_model(Tc) # Run the coalescent model + resolve polytomies
        elif config_dic["resolve_poly"]==True:
            tt.resolve_polytomies() # Resolve multiple mergers
        else:
            pass  # do nothing
    except:
        warn = "Error in the Coalescent model or polytomies resolution."
        _update_session_state(root, state, state_fname=state_fname, warn=warn)


    if config_dic["relax_mu"]==True:  # Relaxed molecular clock 
        _update_session_state(root, state, state_fname=state_fname, err=err)
        # set slack 
        try:
            slack = float(config_dic["slack"])
        except:
            print ("Cannot set mutation rate from the string: %s. Data not understood" % config_dic["mu"])
            slack = None
        if slack <=0 or slack > 1: 
            slack = None

        #set coupling
        try:
            coupling = float(config_dic["coupling"])
        except:
            print ("Cannot set mutation rate from the string: %s. Data not understood" % config_dic["mu"])
            coupling = None
        if coupling <=0 or coupling > 1: 
            coupling = None
        
        if slack is not None and coupling is not None:
            try:
                tt.relaxed_clock(slack, coupling)
            except:
                warn = "Cannot relax molecular clock: exception caught."
                _update_session_state(root, state, state_fname=state_fname, warn=warn)
        else:
            warn = "Cannot relax molecular clock: cannot read slack or coupling variables."
            _update_session_state(root, state, state_fname=state_fname, warn=warn)

    _update_session_state(root, state, state_fname=state_fname, err=err)
    # TODO save results

    _update_session_state(root, state, state_fname=state_fname, err=err)



if __name__=="__main__":
    
    root = "./"
    process(root, config_dic)



