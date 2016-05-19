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
import treetime

dirname = (os.path.dirname(__file__))

def write_json(state_dic, outfile):
    """
    Write pytho dictionary to json file
    """
    with open(outfile, 'w') as of:
            json.dump(state_dic, of)

def create_initial_state(root, config_dic):
    """
    Create the dictionary containing the list of steps and summary of the current
    session state. The computational steps are inferred from the configuration
    dictionary supplied by the end-user.
    """

    # create steps from settings
    state = {
        "state":"Running",
        "warn":[],
        "todo":[],
        "progress":"",
        "done":[],
        "error":""
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
    state["todo"].append("Redirect to the results page")

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

    write_json (state, os.path.join(root, state_fname))
    print (state)
    # TODO save state to json file

def build_tree(root):

    aln_filename = os.path.join(root, "in_aln.fasta")
    tree_filename = os.path.join(root, "in_tree.nwk")
    fast_tree_bin = os.path.join(dirname, "fasttree")
    if not os.path.exists(fast_tree_bin):
        raise (RuntimeError("Cannot find FastTree binary."))

    call = [fast_tree_bin, '-nt','-quiet', aln_filename, ' > ', tree_filename]

    res = os.system(' '.join(call))
    if res != 0:
        raise RuntimeError("FastTree: Exception caught while building the NJ tree")

    return

def save_results(root, tt):
    print (root)
    if not isinstance (tt, treetime.TreeTime):  return

    #  save files
    treetime.treetime_to_json(tt,  os.path.join(root, "out_tree.json"))
    treetime.tips_data_to_json(tt, os.path.join(root, "out_tips.json"))
    treetime.root_lh_to_json(tt,   os.path.join(root, "out_root_lh.json"))
    treetime.root_lh_to_csv(tt,   os.path.join(root, "out_root_lh.csv"))
    # save full alignment
    aln = Align.MultipleSeqAlignment([SeqRecord.SeqRecord (Seq.Seq(''.join(n.sequence)), id=n.name, name=n.name, description="")
        for n in tt.tree.find_clades ()])
    AlignIO.write(aln, os.path.join(root, "out_aln.fasta"), "fasta")

    # save newick tree
    Phylo.write(tt.tree, os.path.join(root, "out_newick_tree.nwk"), "newick")
    #save metadata as csv file
    treetime.save_all_nodes_metadata(tt, os.path.join(root, "out_metadata.csv"))
    #save molecular clock in normal format
    mclock = np.array([(tip.dist2root, tip.numdate_given)
        for tip in tt.tree.get_terminals()
        if hasattr(tip, 'dist2root') and hasattr(tip, 'numdate_given')])
    np.savetxt(os.path.join(root, 'molecular_clock.csv'), mclock,
        delimiter=',',
        header='Distance_to_root,Sampling_date')
    # save GTR in csv file
    treetime.save_gtr_to_file(tt.gtr, os.path.join(root, "out_gtr.txt"))
    # zip all results to one file
    with zipfile.ZipFile(os.path.join(root, 'treetime_results.zip'), 'w') as out_zip:
        out_zip.write(os.path.join(root, 'out_newick_tree.nwk'), arcname='out_newick_tree.nwk')
        out_zip.write(os.path.join(root, 'out_aln.fasta'), arcname='out_aln.fasta')
        out_zip.write(os.path.join(root, 'out_metadata.csv'), arcname='out_metadata.csv')
        out_zip.write(os.path.join(root, 'out_tree.json'), arcname='out_tree.json')
        out_zip.write(os.path.join(root, 'config.json'), arcname='config.json')
        out_zip.write(os.path.join(root, 'molecular_clock.csv'), arcname='molecular_clock.csv')
        out_zip.write(os.path.join(root, 'out_root_lh.csv'), arcname='out_root_lh.csv')
        out_zip.write(os.path.join(root, 'out_gtr.txt'), arcname='out_gtr.txt')


def process(root, config_dic, state_fname="session_state.json"):
    # some global configurations:
    if config_dic["use_mu"]==True:
        try:
            mu = float(config_dic["mu"])
        except Exception, e:
            print ("Cannot set mutation rate from the string. Data not understood" % config_dic["mu"])
            mu = None
    else:
        mu = None

    # specify filenames
    nwk = os.path.join(root, "in_tree.nwk")
    aln = os.path.join(root, "in_aln.fasta")
    meta = os.path.join(root, "in_meta.csv")

    #initialize the session state dictionary from the config_dic:
    state = create_initial_state(root, config_dic)

    if config_dic["do_build_tree"] == True:
        _update_session_state(root, state, state_fname=state_fname)
        try:
            build_tree(root)
        except Exception, e:
            err = "Error occurred while building the tree.  %s" % repr(e)
            _update_session_state(root, state, state_fname=state_fname, err=err)
            return

    _update_session_state(root, state, state_fname=state_fname)

    try:
        gtr = treetime.GTR.standard() # always start from standard J-C model
        tt = treetime.treetime_from_newick(gtr, nwk)  # load tree
        treetime.set_seqs_to_leaves(tt, AlignIO.read(aln, 'fasta'))  # assign sequences
        treetime.read_metadata(tt, meta)  # load meta data
    except Exception, e:
        err = "Cannot initialize TreeTime objects. %s" % repr(e)
        _update_session_state(root, state, state_fname=state_fname, err=err)
        return

    ##  Pipeline
    if config_dic["reuse_branch_len"] == False:
        _update_session_state(root, state, state_fname=state_fname)
        try:
            tt.optimize_seq_and_branch_len(False, True)
        except Exception, e:
            err = "Error in branch len optimization\n%s " % repr(e)
            _update_session_state(root, state, state_fname=state_fname, err=err)
            return


    _update_session_state(root, state, state_fname=state_fname)
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

    except Exception, e:
        err = "Error in TreeTime object initialiation\n%s" %repr(e)
        _update_session_state(root, state, state_fname=state_fname, err=err)
        return

    # run TreeTime-ML optimization
    _update_session_state(root, state, state_fname=state_fname)
    try:
        tt.ml_t()
    except Exception, e:
        err = "Error in TreeTime Max-likelihood optimization\n%s" % repr(e)
        _update_session_state(root, state, state_fname=state_fname, err=err)
        return

    # some post-processing
    _update_session_state(root, state, state_fname=state_fname)
    try:
        if config_dic["coalescent"]==True:
            tt.coalescent_model(Tc) # Run the coalescent model + resolve polytomies
        elif config_dic["resolve_poly"]==True:
            tt.resolve_polytomies() # Resolve multiple mergers
        else:
            pass  # do nothing
    except Exception, e:
        warn = "Error in the Coalescent model run\n%s" % repr(e)
        _update_session_state(root, state, state_fname=state_fname, warn=warn)


    if config_dic["relax_mu"]==True:  # Relaxed molecular clock
        _update_session_state(root, state, state_fname=state_fname)
        # set slack
        try:
            slack = float(config_dic["slack"])
        except Exception, e:
            print ("Cannot set mutation rate from the string\n. Data not understood" % config_dic["mu"])
            slack = None
        if slack <=0 or slack > 1:
            slack = None

        #set coupling
        try:
            coupling = float(config_dic["coupling"])
        except Exception, e:
            print ("Cannot set mutation rate from the string\n. Data not understood" % config_dic["mu"])
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

    _update_session_state(root, state, state_fname=state_fname)
    try:
        save_results(root, tt)
    except Exception, e:
        err = "Error in saving results\n%s" % repr(e)
        _update_session_state(root, state, state_fname=state_fname, err=err)
        return


    _update_session_state(root, state, state_fname=state_fname) # redirection
    _update_session_state(root, state, state_fname=state_fname) # last step -> complete

if __name__=="__main__":

    root = "./"
    process(root, config_dic)



