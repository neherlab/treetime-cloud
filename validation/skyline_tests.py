from external_binaries import FFPOPSIM_SKYLINE_BIN
from utility_functions_simulated_data import _ffpopsim_tree_aln_postprocess, generations_from_ffpopsim_tree
import sys,os, glob

import numpy as np
from treetime import TreeTime
from matplotlib import pyplot as plt
import seaborn as sns

from plot_defaults import *
cols = sns.color_palette(n_colors=3)

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


    call = [FFPOPSIM_SKYLINE_BIN, L, N, SAMPLE_NUM, SAMPLE_FREQ, SAMPLE_VOL, MU, basename, amp, period]
    os.system(' '.join([str(k) for k in call]))

    return basename


def estimate_skyline(base_name, plot=False):
    tree_file = base_name + ".opt.nwk"
    aln_file  = base_name + ".nuc.fasta"
    params = (os.path.split(base_name)[-1]).split('_')
    print(params)
    period = float(params[-2][6:])
    amp = float(params[-3][3:])
    N = float(params[2][1:])

    generations = generations_from_ffpopsim_tree(tree_file)[1]
    T = TreeTime(tree=tree_file, aln=aln_file, dates=generations, gtr="JC69", real_dates=False)

    T.run(Tc="skyline",max_iter=3, long_branch=True, resolve_polytomies=True, infer_gtr=True, root='best') #, fixed_slope=0.0001)
    print(T.gtr)
    skyline = T.merger_model.skyline_inferred()
    skyline_em = T.merger_model.skyline_empirical()


    x = skyline.x
    truePopSize = N*(1.0 + amp*np.cos(2.0*np.pi*x/N/period))
    if plot:
        plt.figure(figsize=onecolumn_figsize)
        plt.plot(x, skyline.y)
        plt.plot(skyline_em.x, skyline_em.y)
        plt.plot(x, truePopSize)

    informative_range = x.searchsorted(np.min([n.numdate for n in T.tree.root]))
    return period, amp, x[informative_range:], skyline.y[informative_range:], truePopSize[informative_range:]

def save_estimate_skyline(period, amp, x, y, trueY, outfile):
    with open(outfile, 'w') as of:
        of.write("#period={}\n".format(period))
        of.write("#amplitude={}\n".format(amp))
        of.write("#x, y, trueY\n")
        for (i,j,k) in zip(x, y, trueY):
            of.write('{},{},{}\n'.format(i,j,k))

def read_estimate_skyline(infile):
    with open(infile, 'r') as inf:
        ss = inf.readlines()

    period = None
    amp = None
    x,y,trueY = [],[],[]

    reading_arrays = False
    for s in ss:
        if s.startswith('#period'):
            period = float(s.split('=')[-1].strip())
        elif s.startswith('#amplitude'):
            amp = float(s.split('=')[-1].strip())
        elif s.startswith('#x, y, trueY'):
            reading_arrays = True
        elif reading_arrays:
            try:
                i,j,k = map(float, s.split(','))
                x.append(i)
                y.append(j)
                trueY.append(k)
            except:
                reading_arrays = False
        else:
            continue
    return period, amp, x, y, trueY

if __name__=="__main__":

    RUN_FFPOPSIM = False
    RUN_TREETIME = False
    SAVE_FIG = True

    sim_dir = './skyline/simulated_data/'
    resdir = './skyline/'
    N=300
    mu = 1e-3
    L=1000
    Nsamples = 40
    DeltaT = N/40
    SampleSize = 20
    periods = [0.5, 2.0]
    bottle_necks = [0.5, 0.8, 0.9]
    if RUN_FFPOPSIM:
        for period in periods:
            for amp in bottle_necks:
                run_ffpopsim_simulation_skyline(L, N, SampleSize, Nsamples, DeltaT, mu, amp, period, sim_dir, 'fluct')


    if RUN_TREETIME:
        fnames = glob.glob(os.path.join(sim_dir,'FF*Mu0.001*fluct.nwk'))
        res = []
        for fname in fnames:
            tmp = estimate_skyline(fname[:-4])
            resname = os.path.join(resdir, os.path.split(fname)[-1].replace('_fluct.nwk', '_res.txt'))
            save_estimate_skyline(tmp[0], tmp[1],tmp[2],tmp[3],tmp[4], resname)
            res.append(tmp)
