#!/usr/bin/env python
import datetime
import utility_functions_flu as flu_utils
import utility_functions_general as gen_utils
import utility_functions_beast as beast_utils

from Bio import AlignIO
import os,sys
import numpy as np

def _run_beast(aln_name, tree_name, known_dates_fraction, out_dir):

        def log_post_process(log_file):

            df = beast_utils.read_beast_log(log_file, np.max(dates.values()))
            if df is None or df.shape[0] < 200:
                print ("Beast log {} is corrupted or BEAST run did not finish".format(log_file))
                return

            inferred_LH = df['likelihood'][-50:].mean()
            inferred_LH_std = df['likelihood'][-50:].std()
            inferred_Tmrca = df['treeModel.rootHeight'][-50:].mean()
            inferred_Tmrca_std = df['treeModel.rootHeight'][-50:].std()
            inferred_Mu = df['clock.rate'][-50:].mean()
            inferred_Mu_std = df['clock.rate'][-50:].std()

            if not os.path.exists(beast_res_file):
                try:
                    with open(beast_res_file, 'w') as of:
                        of.write("#Filename,KnownDatesFraction,LH,LH_std,Tmrca,Tmrca_std,Mu,Mu_std\n")
                except:
                    pass

            with open(beast_res_file, 'a') as of:
                of.write("{},{},{},{},{},{},{},{}\n".format(
                    tree_name,
                    known_dates_fraction,
                    inferred_LH,
                    inferred_LH_std,
                    inferred_Tmrca,
                    inferred_Tmrca_std,
                    inferred_Mu,
                    inferred_Mu_std))

        dates = flu_utils.make_known_dates_dict(aln_name, known_dates_fraction)
        beast_out_dir = os.path.join(out_dir, 'beast_out')
        if not os.path.exists(beast_out_dir):
            try:
                os.makedirs(beast_out_dir)
            except:
                pass

        beast_prefix = os.path.join(beast_out_dir, os.path.split(tree_name)[-1].replace('.nwk', filename_suffix))  # truncate '.nwk'
        flu_utils.run_beast(tree_name, aln_name, dates, beast_prefix,
            log_post_process=log_post_process,
            template_file="./resources/beast/template_bedford_et_al_2015.xml")

if __name__ == "__main__":


    RUN_BEAST = True
    RUN_TREETIME = True

    print sys.argv[0]

    out_dir = sys.argv[1]
    subtree = sys.argv[2]
    known_dates_fraction = float(sys.argv[3])
    filename_suffix= sys.argv[4]


    assert(known_dates_fraction > 0 and known_dates_fraction <= 1.0)

    aln_name = subtree + ".fasta"
    tree_name = subtree + ".nwk"

    if RUN_TREETIME:

        treetime_res_file = os.path.join(out_dir, "treetime_res.csv")

        myTree = flu_utils.create_treetime_with_missing_dates(aln_name, tree_name, known_dates_fraction)
        start = datetime.datetime.now()
        myTree.run(root='best', relaxed_clock=False, max_iter=3, resolve_polytomies=True, do_marginal=False)
        end = datetime.datetime.now()

        if not os.path.exists(treetime_res_file):
            try:
                with open(treetime_res_file, 'w') as of:
                    of.write("#Filename,KnownDatesFraction,Tmrca,Mu,R^2(initial clock),R^2(internal nodes),RunTime(sec)\n")
            except:
                pass

        with open(treetime_res_file, 'a') as of:
            of.write("{},{},{},{},{},{},{}\n".format(
                tree_name,
                str(known_dates_fraction),
                str(myTree.tree.root.numdate),
                str(myTree.date2dist.clock_rate),
                str(myTree.date2dist.r_val),
                str(gen_utils.internal_regress(myTree)),
                str((end-start).total_seconds()) ))

        ##
        ##
        ## precision of the date inference
        ##
        ##
        aln = AlignIO.read(aln_name, 'fasta')
        dates = {k.name: flu_utils.date_from_seq_name(k.name) for k in aln}
        dTs = [(leaf.name, leaf.numdate, dates[leaf.name], leaf.numdate - dates[leaf.name])
                for leaf in myTree.tree.get_terminals() if leaf.numdate_given is None]


        dates_res_file = os.path.join(out_dir, "treetime_dates_res.csv")
        if not os.path.exists(dates_res_file):
            try:
                with open(dates_res_file, 'w') as of:
                    of.write("#LeafName,KnownDatesFraction,Tmrca,LeafDate,LeafDate_real,LeafDateErr(years)\n")
            except:
                pass

        with open(dates_res_file, 'a') as of:
            of.write(
                "\n".join(["{},{},{},{},{},{}".format(
                dT[0],
                str(known_dates_fraction),
                str(myTree.tree.root.numdate),
                str(dT[1]),
                str(dT[2]),
                str(dT[3])) for dT in dTs]) + "\n")


    if RUN_BEAST:
        beast_res_file = os.path.join(out_dir, "beast_res.csv")
        _run_beast(aln_name, tree_name, known_dates_fraction, out_dir)


