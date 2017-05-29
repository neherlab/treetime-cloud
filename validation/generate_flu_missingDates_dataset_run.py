#!/usr/bin/env python
import datetime
import utility_functions_flu as flu_utils
import utility_functions_general as gen_utils

from Bio import AlignIO
import os,sys
if __name__ == "__main__":


    RUN_BEAST = True
    RUN_TREETIME = False

    print sys.argv[0]

    knonw_dates_fraction = float(sys.argv[1])
    out_dir = sys.argv[2]
    aln_name = sys.argv[3]
    tree_name = sys.argv[4]
    res_file = sys.argv[5]
    res_dates_file = sys.argv[6]
    filename_suffix = sys.argv[7]


    assert(knonw_dates_fraction > 0 and knonw_dates_fraction <= 1.0)

    if RUN_TREETIME:
        myTree = flu_utils.create_treetime_with_missing_dates(aln_name, tree_name, knonw_dates_fraction)
        start = datetime.datetime.now()
        myTree.run(root='best', relaxed_clock=False, max_iter=3, resolve_polytomies=True, do_marginal=False)
        end = datetime.datetime.now()
        with open(res_file, 'a') as of:
            of.write("{},{},{},{},{},{},{}\n".format(
                tree_name,
                str(knonw_dates_fraction),
                str(myTree.tree.root.numdate),
                str(myTree.date2dist.slope),
                str(myTree.date2dist.r_val),
                str(gen_utils.internal_regress(myTree)),
                str((end-start).total_seconds()) ))

        # precision of the date inference
        aln = AlignIO.read(aln_name, 'fasta')
        dates = {k.name: flu_utils.date_from_seq_name(k.name) for k in aln}
        dTs = [(leaf.name, leaf.numdate, dates[leaf.name], leaf.numdate - dates[leaf.name])
                for leaf in myTree.tree.get_terminals() if leaf.numdate_given is None]

        with open(res_dates_file, 'a') as of:
            of.write(
                "\n".join(["{},{},{},{},{},{}".format(
                dT[0],
                str(knonw_dates_fraction),
                str(myTree.tree.root.numdate),
                str(dT[1]),
                str(dT[2]),
                str(dT[3])) for dT in dTs]) + "\n")

    if RUN_BEAST:
        dates = flu_utils.make_known_dates_dict(aln_name, knonw_dates_fraction)
        beast_out_dir = os.path.join(out_dir, 'beast_out')
        if not os.path.exists(beast_out_dir):
            try:
                os.makedirs(beast_out_dir)
            except:
                pass
        beast_prefix = os.path.join(beast_out_dir, os.path.split(tree_name)[-1].replace('.nwk', filename_suffix))  # truncate '.nwk'
        flu_utils.run_beast(tree_name, aln_name, dates, beast_prefix, template_file="./resources/beast/template_bedford_et_al_2015.xml")


