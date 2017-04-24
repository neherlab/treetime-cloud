from external_binaries import FFPOPSIM_SKYLINE_BIN

from generate_dataset import _ffpopsim_tree_aln_postprocess,
import sys,os, glob

import numpy as np
from treetime import TreeTime
from treetime import io
from generate_dataset import dates_from_ffpopsim_tree
from matplotlib import pyplot as plt
import seaborn as sns

cols = sns.color_palette(n_colors=6)

def run_ffpopsim_simulation_skyline(L, N, SAMPLE_VOL, SAMPLE_NUM, SAMPLE_FREQ, MU, amplitude,
                            period, res_dir, res_suffix, failed=None, **kwargs):
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
    basename = _run_ffpopsim_skyline(L=L, N=N,
                    SAMPLE_NUM=SAMPLE_NUM,
                    SAMPLE_FREQ=SAMPLE_FREQ,
                    SAMPLE_VOL=SAMPLE_VOL,
                    MU=MU, amp=amplitude, period=period,
                    res_dir=res_dir, res_suffix=res_suffix)

    # pot-process the results
    if 'optimize_branch_len' in kwargs:
        optimize_branch_len = kwargs['optimize_branch_len']
    else:
        optimize_branch_len = True
    _ffpopsim_tree_aln_postprocess(basename, optimize_branch_len=optimize_branch_len)
    print ("Done clusterSingleFunc")
    return basename

def _run_ffpopsim_skyline(L=100, N=100, SAMPLE_NUM=10, SAMPLE_FREQ=5, SAMPLE_VOL=15, MU=5e-5,
                          amp=0.9, period = 1.0, res_dir="./", res_suffix=""):
    """
    Simple wrapper function to call FFpopSim binary in a separate subprocess
    """

    basename = "FFpopSim_L{}_N{}_Ns{}_Ts{}_Nv{}_Mu{}_Amp{}_Tfluct{}".format(str(L), str(N),
                str(SAMPLE_NUM), str(SAMPLE_FREQ), str(SAMPLE_VOL), str(MU), str(amp), str(period))

    basename = os.path.join(res_dir, basename)
    if res_suffix != "":
        basename = basename + "_" + res_suffix


    call = [FFPOPSIM_BIN, L, N, SAMPLE_NUM, SAMPLE_FREQ, SAMPLE_VOL, MU, basename, amp, period]
    os.system(' '.join([str(k) for k in call]))

    return basename


def estimate_skyline(base_name, plot=False):
    tree_file = "out/%s.opt.nwk"%base_name
    aln_file  = "out/%s.nuc.fasta"%base_name
    params = base_name.split('_')
    print(params)
    period = float(params[-2][6:])
    amp = float(params[-3][3:])
    N = float(params[2][1:])

    T = TreeTime(tree=tree_file, aln=aln_file, dates = dates_from_ffpopsim_tree(tree_file)[1], gtr="JC69", real_dates=False)

    T.run(Tc="skyline",max_iter=3, long_branch=True, resolve_polytomies=True, infer_gtr=True, root='best') #, fixed_slope=0.0001)
    print(T.gtr)
    skyline = T.merger_model.skyline_inferred()
    skyline_em = T.merger_model.skyline_empirical()


    x = skyline.x
    truePopSize = N*(1.0 + amp*np.cos(2.0*np.pi*x/N/period))
    if plot:
        plt.figure()
        plt.plot(x, skyline.y)
        plt.plot(skyline_em.x, skyline_em.y)
        plt.plot(x, truePopSize)

    informative_range = x.searchsorted(np.min([n.numdate for n in T.tree.root]))
    return period, amp, x[informative_range:], skyline.y[informative_range:], truePopSize[informative_range:]


if __name__=="__main__":

    res_dir = './skyline/simulated_data/'
    N=300
    mu = 1e-3
    L=1000
    Nsamples = 40
    DeltaT = N/40
    SampleSize = 20
    periods = [0.5, 1.0, 2.0]
    for period in periods:
        for amp in [0.5, 0.8, 0.9]:
            run_ffpopsim_simulation_skyline(L, N, SampleSize, Nsamples, DeltaT, mu, amp, period, res_dir, 'fluct')

    fnames = glob.glob(os.path.join(res_dir,'FF*Mu0.001*fluct.nwk'))
    res = []
    for fname in fnames:
        tmp = estimate_skyline('.'.join(fname.split('/')[-1].split('.')[:-1]))
        res.append(tmp)

    fig, axs = plt.subplots(3,1)
    for ri, (p, amp, x, s, t) in enumerate(res):
        ax = axs[periods.index(p)]
        print(p,amp, np.corrcoef(s,t)[0,1])
        ax.plot(x, s, c=cols[ri%len(cols)], ls='--')
        ax.plot(x, t, c=cols[ri%len(cols)], ls='-')


