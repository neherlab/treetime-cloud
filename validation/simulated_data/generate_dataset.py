"""
This script defines functions to produce dataset and to run LSD/TreeTime
reconstruction on the given dataset.
The results of the LSD/TReeTime are saved into given files for further analysis
and visualizations.

The functions in this script depend on the external tools.
Define the paths to the FFPopSim, FastTree, LSD binaries in order to run the full set
of the simulations

The examples of the dataset generation can be found in a separate file.
"""

from __future__ import print_function, division

from Bio import Phylo, Seq, SeqRecord, Align, AlignIO
import os, sys
import subprocess
import numpy as np
import treetime
from scipy.stats import linregress

FFPOPSIM_BIN = "/ebio/ag-neher/share/users/psagulenko/Programming/FFPopSim/FFpopSim_exe/FFpopSim_exe" #"/home/psagulenko/PHD/prog/FFpopSim-QT/build-FFpopSim-Desktop-Debug/FFpopSim"
FAST_TREE_BIN = "/ebio/ag-neher/share/users/psagulenko/programs/fast_tree/fasttree" #"/home/psagulenko/PHD/programs/fast_tree/fasttree"
LSD_BIN = "/ebio/ag-neher/share/users/psagulenko/programs/LSD/lsd"

def run_ffpopsim_simulation(L, N, SAMPLE_VOL, SAMPLE_NUM, SAMPLE_FREQ, MU, res_dir, res_suffix, failed=None, **kwargs):
    """
    Run simulation with FFPopSim package and perform the data preprocessing.
    The FFpopSim produces phylogenetic tree and alignment in the binary (0-1) form.
    The tree branch lengths are in units of time, expressed in generations.
    Returns:
     - basename: base name of the files, where the results are stored. The file
     suffixes are added for each file type separately.
    """

    sys.stdout.write("Importing modules...")

    # check the output location
    if not os.path.exists(res_dir) or not os.path.isdir(res_dir):
        os.makedirs(res_dir)

    # run ffpopsim
    sys.stdout.write("Running FFpopSim...")
    basename = _run_ffpopsim(L=L, N=N,
                    SAMPLE_NUM=SAMPLE_NUM,
                    SAMPLE_FREQ=SAMPLE_FREQ,
                    SAMPLE_VOL=SAMPLE_VOL,
                    MU=MU,
                    res_dir=res_dir, res_suffix=res_suffix)

    # pot-process the results
    if 'optimize_branch_len' in kwargs:
        optimize_branch_len = kwargs['optimize_branch_len']
    else:
        optimize_branch_len = True
    _ffpopsim_tree_aln_postprocess(basename, optimize_branch_len=optimize_branch_len)
    print ("Done clusterSingleFunc")
    return basename

def _run_ffpopsim(L=100, N=100, SAMPLE_NUM=10, SAMPLE_FREQ=5, SAMPLE_VOL=15, MU=5e-5, res_dir="./", res_suffix=""):
    """
    Simple wrapper function to call FFpopSim binary in a separate subprocess
    """

    basename = "FFpopSim_L{}_N{}_Ns{}_Ts{}_Nv{}_Mu{}".format(str(L), str(N),
                str(SAMPLE_NUM), str(SAMPLE_FREQ), str(SAMPLE_VOL), str(MU))

    basename = os.path.join(res_dir, basename)
    if res_suffix != "":
        basename = basename + "_" + res_suffix


    call = [FFPOPSIM_BIN, L, N, SAMPLE_NUM, SAMPLE_FREQ, SAMPLE_VOL, MU, basename]
    os.system(' '.join([str(k) for k in call]))

    return basename

def _ffpopsim_tree_aln_postprocess(basename, optimize_branch_len=False, prefix='Node/'):

    """
    Given the raw data produced ini the FFPopSim simulation, perform the preliminary
    data processing. This includes converting the alignment in nucleotide notation, where
    0 is converted to 'A', 1 is converted to 'C'. The tree branch lengths are been
    optimized. The tree is checked to have no multiple mergers (the multiple mergers are
    resolved randomly). The nodes of tree are named in a unified way and if there are
    non-unique names found, they resolved by adding additional suffixes
    """

    def ffpopsim_aln_to_nuc(basename):

        with open(basename + ".bin.fasta", 'r') as inf:
            ss = inf.readlines()

        with open(basename + ".nuc.fasta", 'w') as of:
            for s in ss:
                if s.startswith(">"):
                    of.write(s)
                else:
                    s = s.replace('0', "A").replace('1', "C")
                    of.write(s)

    def generation_from_node_name(name):
        try:
            return int(name.split("_")[1])
        except:
            return -1


    from collections import Counter

    t = Phylo.read(basename + ".nwk", "newick")
    t = _remove_polytomies(t)
    cs = [k for k in t.find_clades() if len(k.clades)==1]
    for clade in cs:
        t.collapse(clade)
    t.root.branch_length = 1e-5
    t.ladderize()
    # remove duplicate names:
    t_counter = Counter(t.find_clades())
    for clade in t.find_clades():
        if clade.name is None:
            continue
        if t_counter[clade] > 1 :
            name_suffix += t_counter[clade]
            clade.name += "/"+str(name_suffix)
            t_counter[clade] -= 1

        if clade.name is not None and not clade.name.startswith(prefix):
            clade.name = prefix + clade.name

    max_generation = np.max([generation_from_node_name(k.name) for k in t.find_clades()])
    min_generation = generation_from_node_name(t.root.name)
    assert(max_generation > min_generation and min_generation > 0)

    for clade in t.find_clades():
        node_gen = generation_from_node_name(clade.name)
        if clade.name is not None and not '_DATE_' in clade.name and node_gen != -1:

            clade.name += "_DATE_" + str(2016.5 - (max_generation - node_gen))


    Phylo.write(t, basename + ".nwk", "newick")

    # prepare alignment
    ffpopsim_aln_to_nuc(basename)
    aln = AlignIO.read(basename + ".nuc.fasta", "fasta")
    names = [k.id for k in aln]
    aln_counter = Counter(names)
    for k in aln:
        key = k.id
        # remove duplicate names
        if (aln_counter[key] > 1):
            k.id = k.id + "/" + str(aln_counter[k.id])
            k.name = k.id
            aln_counter[key] -= 1
        # add name prefix
        if not k.id.startswith(prefix):
            k.id = prefix + k.id
            k.name = k.id

        node_gen = generation_from_node_name(k.id)
        if not "_DATE_" in k.id and node_gen != -1:
            k.id += "_DATE_" + str(2016.5 - (max_generation - node_gen))
            k.name = k.id

    AlignIO.write(aln, basename + ".nuc.fasta", "fasta")

    if optimize_branch_len:
        import treetime
        gtr = treetime.GTR.standard()
        tanc = treetime.TreeAnc(aln=aln,tree=t,gtr=gtr)
        tanc.optimize_seq_and_branch_len(reuse_branch_len=False,prune_short=False,infer_gtr=False)
        Phylo.write(tanc.tree, basename+".opt.nwk", "newick")

def _remove_polytomies(tree):
    """
    Scan tree, and remove the polytomies (if any) by random merging
    Args:
     - tree: Phylogenetic tree as Bio.Phylo object
    Returns: tree without polytomies
    """
    for clade in tree.find_clades():
        #if not hasattr(clade, "name") or clade.name is None:
        #    clade.name = "None"
        if len(clade.clades) < 3:
            continue
        clades = clade.clades
        while len(clades) > 2:
            c1 = clades.pop()
            c2 = clades.pop()
            new_node = Phylo.BaseTree.Clade()
            new_node.name=None
            new_node.branch_length = 1e-6
            new_node.clades = [c1,c2]
            clades.append(new_node)
        clade.clades = clades

    return tree

def reconstruct_fasttree(basename, optimize_branch_len=False):
    """
    Run the fast tree reconstruction given the alignment produced by FFPopSim
    simulations
    """
    def fasttree_post_process(aln, basename):
        #import treetime
        ffpopsim_treefile = basename + ".nwk"
        treefile = basename + ".ft.nwk"
        alnfile = basename + ".nuc.fasta"
        tree = Phylo.read(treefile, 'newick')
        tree = _remove_polytomies(tree)
        tree.ladderize()
        Tmrca,dates = _dates_from_ffpopsim_tree(ffpopsim_treefile)
        tree.root.name = "FFPOPsim_Tmrca_DATE_" + str(Tmrca)
        return tree
        ## optimize branch lengths (mainly, to remove the fasttree zero-lengths artefacts)
        #if optimize_branch_len:
        #    gtr = treetime.GTR.standard()
        #    tanc = treetime.TreeAnc(tree=tree, aln=aln, gtr=gtr)
        #    tanc.optimize_seq_and_branch_len(reuse_branch_len=True,prune_short=False,infer_gtr=False)
        #    return tanc.tree
        #else:
        #    return tree
    fasta = basename + ".nuc.fasta"
    outfile = fasta.replace('.nuc.fasta', ".ft.nwk")
    # call FastTree to reconstruct newick tree
    call = [FAST_TREE_BIN, "-nt", fasta, ">", outfile]
    os.system(' '.join([str(k) for k in call]))
    #subprocess.call(FAST_TREE_BIN +  " -nt " + fasta + " > " + outfile, shell=False)
    tree = fasttree_post_process(fasta, basename)
    os.remove(outfile)
    Phylo.write(tree, outfile, 'newick')

def _dates_from_ffpopsim_tree(t):
    """
    Args:
     - t: tree filename or BioPyhton tree object
    """
    if isinstance(t, str):
        t = Phylo.read(t, 'newick')
    try:
        dates = {}
        for clade in t.get_terminals():
            if "_DATE_"not in clade.name:
                continue
            dates[clade.name] = float(clade.name.split("_")[-1])
        Tmrca = float (t.root.name.split("_")[-1])
        return Tmrca, dates
    except:
        return 2016.5, {}

def run_lsd(treefile, datesfile, outfile, res_file):
    """
    Infer the dates of the internal nodes using the LSD package.
    Append results to the given file.
    Args:

     - treefile: the name for the input tree in newick format.

     - datesfile:  Filename, where to store the dates in the LSD format. If the
     file does not exist, it will be created automatically.

     - outfile: prefix for the filename, where LSD will store the results (different
     trees, output logs, etc.).

     - res_file: file where to write the formatted result string. This file is
     then parsed by the processing scripts to produce comparison with TreeTime
    """
    outfile = outfile.replace(".nwk", ".res.txt")
    Tmrca, dates = _create_date_file_from_ffpopsim_tree(treefile, datesfile)
    # call LSD binary
    call = [LSD_BIN, '-i', treefile, '-d', datesfile, '-o', outfile,
            '-r']
    subprocess.call(call)
    # parse LSD results
    try:
        with open (outfile, 'r') as inf:
            ss = inf.readlines()
    except:
        with open (outfile + 'q', 'r') as inf:
            ss = inf.readlines()
    mystr = [i for i in ss if 'tMRCA' in i][0]
    tmrca = (mystr.split(" ")[5])
    mu = (mystr.split(" ")[3])
    objective = (mystr.split(" ")[7])
    with open(res_file, "a") as of:
        of.write(",".join([treefile, str(Tmrca), tmrca, mu, objective]))

def _create_date_file_from_ffpopsim_tree(treefile, datesfile):
    """
    Read dates of the terminal nodes in the simulated tree, and save them into the
    dates file in LSD format
    """
    Tmrca, dates = _dates_from_ffpopsim_tree(treefile)
    with open(datesfile, 'w') as df:
        df.write(str(len(dates)) + "\n")
        df.write("\n".join([str(k) + "\t" + str(dates[k]) for k in dates]))
    return Tmrca, dates

def run_treetime(basename, outfile, fasttree=False, failed=None):
    """
    Infer the dates of the internal nodes using the TreeTime package.
    Append results to the given file.
    """
    def internal_regress(myTree):
        resarr = []
        for node in myTree.tree.get_nonterminals():
            try:
                resarr.append((node.numdate, node.dist2root))
            except:
                continue
        resarr = np.array(resarr)
        if resarr.shape[0] == 0:
            return 0.
        else:
            return linregress(resarr[:, 0], resarr[:, 1]).rvalue**2

    try:
        if fasttree:
            treefile = basename + ".ft.nwk"
        else:
            treefile = basename + ".opt.nwk"
        aln = basename+'.nuc.fasta'
        Tmrca, dates = _dates_from_ffpopsim_tree(Phylo.read(treefile, "newick"))
        myTree = treetime.TreeTime(gtr='Jukes-Cantor', tree = treefile,
            aln = aln, verbose = 4, dates = dates, debug=False)
        if not fasttree:
            myTree.infer_ancestral_sequences(method='fitch')
            myTree.optimize_branch_len()
            root = None
        else:
            root = 'best'
        myTree.run(root=root, relaxed_clock=False, max_iter=1,resolve_polytomies=False, do_marginal=False)
        with open(outfile, 'a') as of:
            of.write("{},{},{},{},{},{}\n".format(
                basename,
                str(Tmrca),
                str(myTree.tree.root.numdate),
                str(myTree.date2dist.slope),
                str(myTree.date2dist.r_val),
                str(internal_regress(myTree)) ))
        return myTree
    except:
        if failed is not None:
            failed.append(basename)

def evolve_seq(treefile, basename, mu=0.0001, L=1000, mygtr = treetime.GTR.standard()):
    """
    Generate a random sequence of a given length, and evolve it on the tree

    Args:
     - treefile: filename for the tree, on which a sequence should be evolved.
     - basename: filename prefix to save alignments.
     - mu: mutation rate. The units of the mutation rate should be consistent with
     the tree branch length
     - L: sequence length.
     - mygtr: GTR model for sequence evolution

    """
    from treetime import TreeAnc, GTR
    from treetime import seq_utils
    from Bio import Phylo, AlignIO
    from StringIO import StringIO
    import numpy as np
    from itertools import izip

    mygtr.mu = mu
    tree = Phylo.read(treefile, 'newick')
    tree.root.ref_seq = np.random.choice(mygtr.alphabet, p=mygtr.Pi, size=L)
    print ("Started sequence evolution...")
    mu_real = 0.0
    n_branches = 0
    #print ("Root sequence: " + ''.join(tree.root.ref_seq))
    for node in tree.find_clades():
        for c in node.clades:
            c.up = node
        if hasattr(node, 'ref_seq'):
            continue
        t = node.branch_length
        p = mygtr.propagate_profile( seq_utils.seq2prof(node.up.ref_seq, mygtr.profile_map), t)
        # normalie profile
        p=(p.T/p.sum(axis=1)).T
        # sample mutations randomly
        ref_seq_idxs = np.array([int(np.random.choice(np.arange(p.shape[1]), p=p[k])) for k in np.arange(p.shape[0])])
        node.ref_seq = np.array([mygtr.alphabet[k] for k in ref_seq_idxs])
        node.ref_mutations = [(anc, pos, der) for pos, (anc, der) in
                            enumerate(izip(node.up.ref_seq, node.ref_seq)) if anc!=der]
        #print (node.name, len(node.ref_mutations))
        mu_real += 1.0 * (node.ref_seq != node.up.ref_seq).sum() / L
        n_branches += t
    mu_real /= n_branches
    print ("Mutation rate is {}".format(mu_real))
    records = [Align.SeqRecord(Align.Seq("".join(k.ref_seq)), id=k.name, name=k.name)
        for k in tree.get_terminals()]
    full_records = [Align.SeqRecord(Align.Seq("".join(k.ref_seq)), id=k.name, name=k.name)
        for k in tree.get_terminals()]
    #import ipdb; ipdb.set_trace()
    aln = Align.MultipleSeqAlignment(records)
    full_aln = Align.MultipleSeqAlignment(full_records)
    print ("Sequence evolution done...")

    # save results
    AlignIO.write(aln, basename+'.aln.ev.fasta', 'fasta')
    AlignIO.write(full_aln, basename+'.aln.ev_full.fasta', 'fasta')

    return aln, full_aln, mu_real

if __name__=="__main__":

    pass
    #df_tt = get_datatreetime_frame('./accuracy_4_treetime_res.csv')
    #df_tt_ft = get_datatreetime_frame('./accuracy_4_treetime_res_ftree.csv')
    #df_lsd = get_LSD_frame('./accuracy_4_lsd_res.csv')
    #df_lsd_ft = get_LSD_frame('./accuracy_4_lsd_res_ftree.csv')
#
#    #plot_dataframe(df_tt, tt_or_lsd="TreeTime", TN_threshold=10, FAST_TREE=False,CLEAR_FIGURE=False)
#    #plot_dataframe(df_tt_ft, tt_or_lsd="TreeTime", TN_threshold=10, FAST_TREE=True,CLEAR_FIGURE=False)
#    #plot_dataframe(df_lsd, tt_or_lsd="LSD", TN_threshold=10, FAST_TREE=False,CLEAR_FIGURE=False)
#    #plot_dataframe(df_lsd_ft, tt_or_lsd="LSD", TN_threshold=10, FAST_TREE=True,CLEAR_FIGURE=False)
#    ## df_tt = get_datatreetime_frame("./accuracy_3_res.csv")
#    ## df_lsd = get_LSD_frame("./accuracy_3_lsd_res.csv")
#
#    #plt.figure(4)
#    #plt.xscale('log')
#    #plt.figure(5)
    #plt.xscale('log')


















