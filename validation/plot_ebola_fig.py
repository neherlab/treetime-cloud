'''
this script illustrates the use of treetime to analyze viral sequences.
As an example, we use Ebola virus sequences from the 2014-2015 outbreak in
West Africa. This example contains more than 300 sequences, it will take
a few minutes to run.
'''

from __future__ import print_function, division
from treetime import TreeTime
import numpy as np
from scipy import optimize as sciopt

fs=16


def load_case_numbers(res=5):
    w = np.ones(res, dtype=float)/res
    import pandas as pd
    from datetime import datetime
    data = pd.read_csv('./resources/ebola/ebola_cumulative_reported_cases.csv')
    total = np.diff(np.sum([data.loc[:,x] for x in data.columns if 'ases' in x], axis=0))
    dates = [datetime.strptime(x, "%m/%d/%Y") for x in data.loc[:,"WHO report date"]]
    numdates = [x.year + x.timetuple().tm_yday/365.25 for x in dates]
    return np.convolve(w, numdates[:-1], mode='valid'), np.convolve(w, total, mode='valid')

if __name__ == '__main__':

    # load data and parse dates
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import seaborn as sns
    sns.set_style('whitegrid')
    from Bio import Phylo
    plt.ion()
    base_name = './resources/ebola/ebola'
    import datetime
    #from treetime.utils import numeric_date
    with open(base_name+'.csv') as date_file:
        dates = {}
        for line in date_file:
            if line[0]=='#':
                continue
            try:
                name, date = line.strip().split(',')
                dates[name] = float(date)
            except:
                continue

    # instantiate treetime
    ebola = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                        aln = base_name+'.fasta', verbose = 4, dates = dates)

    # infer an ebola time tree while rerooting and resolving polytomies
    ebola.run(root='best', relaxed_clock=False, max_iter=2,
              resolve_polytomies=True, Tc='skyline', do_marginal="assign")
    # plot Skyline
    skyline, confidence = ebola.merger_model.skyline_inferred(gen=50, confidence=2.0)
    skyline_empirical = ebola.merger_model.skyline_empirical(gen=50)

    # scatter root to tip divergence vs sampling date
    ebola.plot_root_to_tip(add_internal=True)
    t=np.array([2014,2016])
    plt.plot(t, t*ebola.date2dist.clock_rate + ebola.date2dist.intercept,
             label="y = %1.5f t%1.3f"%(ebola.date2dist.clock_rate, ebola.date2dist.intercept))
    plt.legend(loc=2)

    # rescale branch length to years and plot in axis 0
    from treetime import plot_vs_years
    fig, axs = plt.subplots(2,1, sharex=True, figsize=(10,12))
    plot_vs_years(ebola, years=1, ax=axs[0], confidence=(0.05,0.95), label_func = lambda x:"")
    axs[0].set_xlim(0, 2.5)
    axs[0].tick_params(labelsize=fs*0.8)

    # reset branch length to time (in substitution rate units)
    for n in ebola.tree.find_clades():
        if n.up:
            n.branch_length=n.clock_length


    axs[1].fill_between(skyline.x-ebola.tree.root.numdate, confidence[0], confidence[1], color=(0.8, 0.8, 0.8))
    axs[1].plot(skyline.x-ebola.tree.root.numdate, skyline.y, label='coalescent population size estimate')
    dates, cases = load_case_numbers(res=5)
    axs[1].plot(dates-ebola.tree.root.numdate, cases, label='WHO report')
    axs[1].plot([0,2.5], [1,1])
    # axs[1].plot(skyline_empirical.x, skyline_empirical.y, label='empirical skyline')
    axs[1].set_ylim([0.5, 1200])
    axs[1].set_yscale('log')
    axs[1].set_xlabel('date', fontsize=fs)
    axs[1].tick_params(labelsize=fs*0.8)
    axs[1].legend(fontsize=fs)
    plt.tight_layout()

    plt.savefig('figs/ebola.pdf')
    plt.savefig('figs/ebola.png')
    plt.savefig('figs/ebola.svg')




