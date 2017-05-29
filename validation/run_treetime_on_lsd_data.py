from Bio import AlignIO, Phylo
from treetime_repo import TreeTime
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np

def load_dates(data_set):
    dates = {}
    with open('LSD_validation_data/true_trees/%s.date'%data_set) as ifile:
        for line in ifile:
            entries = line.strip().split()
            if len(entries)==2:
                dates[entries[0]]=float(entries[1])
    return dates

def plot_stretch(tt):
    fig = plt.figure()
    ax = plt.subplot(111)
    vmin, vmax = 0.5, 1.5 # color branches according to the rate deviation
    for n in tt.tree.find_clades():
        if n.up:
            n.color = [int(x*255) for x in cm.cool((min(max(vmin, n.branch_length_interpolator.gamma),vmax)-vmin)/(vmax-vmin))[:3]]
        else:
            n.color = [200,200,200]

    ax.set_title("relaxed clock")
    Phylo.draw(tt.tree, axes=ax, show_confidence=False,
                                label_func = lambda x:'')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            description="Run treetime on simulated data used to test LSD in To et al")
    parser.add_argument('--tree', required = True, type = str,  help ="Which tree to use, one of D750_3_25, D750_11_10, D995_11_10, D995_3_25")
    parser.add_argument('--tree_method', type = str, default="true" help ="Which reconstruction method to use. one of true or Phyml")
    parser.add_argument('--tree_rooting', type = str, default="unrooted" help ="Use rooted trees or not")
    parser.add_argument('--clock', type = str, default='strict',
                        help ="strict or relaxed molecular clock")
    parser.add_argument('--relax',nargs='*', default = False,
                        help='use an autocorrelated molecular clock. Prior strength and coupling of parent '
                             'and offspring rates can be specified e.g. as --relax 1.0 0.5')

    args = parser.parse()
    data_set = args.tree
    rate_type = args.clock
    # specify relaxed clock if desired -> make dict to pass to treetime
    if args.relax:
        rc = {'slack':float(args.relax[0]), 'coupling':float(args.relax[1])}
    else:
        rc = False

    # convert input arguments to the method of tree reconstruction which differs by type
    # currently only Phyml and true trees included
    if args.tree_method=='true':
        tree_set = ['true', 'out']
    elif args.tree_method=='Phyml':
        tree_set = ['Phyml', args.tree_rooting]

    # true rate
    rate = 0.006

    # load simulated data
    trees = list(Phylo.parse('LSD_validation_data/%s/%s/%s_%s.tree'%(tree_set[0], rate_type, data_set,tree_set[1]), 'newick'))
    alns = list(AlignIO.parse('LSD_validation_data/Alignments/%s/%s_out.phy'%(rate_type, data_set), 'phylip'))
    dates = load_dates(data_set)

    # loop over 100 trees and collect treetime results
    res = []
    for ti in range(len(trees)):
        T = trees[ti]
        out_groups = [n for n in T.get_terminals() if n.name=='out']
        if len(out_groups):
            T.prune(out_groups[0])
        tt = TreeTime(tree=T, aln=alns[ti], dates=dates, gtr='JC69')
        tt.run(root='best', infer_gtr=True, max_iter=1, n_iqd=4,
               relaxed_clock=rc, use_input_branch_length=True)
        div = [n.branch_length for n in tt.tree.find_clades() if n.up]
        W = tt.gtr.W
        # this stores the clock rate, the root date, the average branch length
        # GTR.Pi and the transition/transversion rates
        res.append([tt.date2dist.slope, tt.tree.root.numdate, np.mean(div)]
                     + list(tt.gtr.Pi[:4]) + [W[0,2], W[1,3]] +
                    [(W[0,1]+W[0,3]+W[1,2]+W[2,3])/4.0])


    res = np.array(res)
    rate_bias = np.mean(res[:,0]-rate)
    rate_error = np.sqrt(np.mean((res[:,0]-rate)**2))
    print("rate bias: %0.5f, error: %1.5f"%(rate_bias/rate, rate_error/rate))

    TMRCA_bias = np.mean(res[:,1])
    TMRCA_error = np.sqrt(np.mean((res[:,1])**2))
    print("TMRCA bias: %0.5f, error: %1.5f"%(TMRCA_bias, TMRCA_error))

    # write to file
    np.savetxt('LSD_validation_data/TT_%s_%s_%s.txt'%(data_set, rate_type,
                                    ('rc' if rc else 'sc')), res)
