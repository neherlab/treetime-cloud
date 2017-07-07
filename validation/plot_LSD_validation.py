import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

from plot_defaults import *
save_fig = True

parser = argparse.ArgumentParser(
        description="Run treetime on simulated data used to test LSD in To et al")
parser.add_argument('--tree', required = True, type = str,  help ="Which tree to use, one of 750_3_25, 750_11_10, 995_11_10, 995_3_25")

args = parser.parse_args()

data_set = args.tree
data = pd.read_csv('LSD_validation_data/Simulation_results_%s.csv'%data_set)


rate_model = {}
method_map = {}
topo_map = {}
last_method = ''
last_rate = ''
for i in range(51):
    if type(data.iloc[0,i])==str:
        last_rate = data.iloc[0,i].split()[0].lower()
    if type(data.iloc[1,i])==str:
        last_method = data.iloc[1,i]
    method_map[i]=last_method
    topo_map[i]=data.iloc[2,i]
    rate_model[i]=last_rate


rates = {}
TMRCAs = {}

for i in range(3,51):
    rates[(rate_model[i], method_map[i], topo_map[i])] = np.array(data.iloc[3:103,i].astype(float))
    TMRCAs[(rate_model[i], method_map[i], topo_map[i])] = np.array(data.iloc[103:203,i].astype(float))

rate=0.006
for model in rates:
    print(str(model)+"rate error: %05f, bias: %05f"%(np.sqrt(np.mean((rates[model]-rate)**2))/rate, np.mean(rates[model]-rate)/rate))

print("****\n\n****")
for model in TMRCAs:
    print(str(model)+"TMRCA error: %05f, bias: %05f"%(np.sqrt(np.mean((TMRCAs[model])**2)), np.mean(TMRCAs[model])))

TT_data={}
for rate_type in ["strict", "relaxed"]:
    rc="sc"
    label = "D%s_%s_%s"%(data_set, rate_type, rc)
    TT_data[label] = np.loadtxt('LSD_validation_data/TT_%s.txt'%label)
    print(label, TT_data[label].mean(axis=0))
    bl_factor = 1.0 #np.exp(4.*TT_data[label].mean(axis=0)[2])
    tree_type = "PhyML"
    methods = set([x for x in method_map.values() if len(x) and x[-1]!='*'])

    tmp = {'TT':(TT_data[label][:,0]*bl_factor-rate)/rate,
            'BSMC':(rates[(rate_type, "BSMC", "True\ntopology")]-rate)/rate}
    # if rate_type=='relaxed':
    #     tmp['BRMC'] = (rates[(rate_type, "BRMC", "True\ntopology")]-rate)/rate
    tmp.update({m:(rates[(rate_type, m, tree_type)]-rate)/rate for m in methods if (rate_type, m, tree_type) in rates})
    df = pd.DataFrame(tmp)

    fig = plt.figure(figsize=onecolumn_figsize)
    axes = fig.add_subplot(111)
    axes.hlines(0, -1, len(df), lw=3)
    sns.violinplot(df)
    #axes.set_title('Evolution model: %s'%rate_type)
    plt.ylabel('relative clock rate error, $[\Delta \mu/\mu]$', fontsize=label_fs)
    for tick_label in axes.get_xticklabels():
            tick_label.set_fontsize(label_fs)
    for tick_label in axes.get_yticklabels():
            tick_label.set_fontsize(tick_fs)

    if save_fig:
        for fmt in ['pdf', 'png', 'svg']:
            fig.savefig("./figs/LSD_data_%s_%s_Mu.%s"%(data_set, rate_type, fmt))

    tmp = {'TT':TT_data[label][:,1],
            'BSMC':TMRCAs[(rate_type, "BSMC", "True\ntopology")]}
    if rate_type=='relaxed':
        tmp['BRMC'] = (rates[(rate_type, "BRMC", "True\ntopology")]-rate)/rate
    tmp.update({m:TMRCAs[(rate_type, m, tree_type)] for m in methods if (rate_type, m, tree_type) in TMRCAs})
    df = pd.DataFrame(tmp)

    fig = plt.figure(figsize=onecolumn_figsize)
    axes = fig.add_subplot(111)
    # axes.set_title('Evolution model: %s'%rate_type)
    axes.hlines(0, -1, len(df), lw=3)
    sns.violinplot(df)
    axes.set_ylabel('TMRCA error [$\mathrm{Years}$]', fontsize=label_fs)
    for tick_label in axes.get_xticklabels():
            tick_label.set_fontsize(label_fs)
    for tick_label in axes.get_yticklabels():
            tick_label.set_fontsize(tick_fs)

    if save_fig:
        for fmt in ['pdf', 'png', 'svg']:
            fig.savefig("./figs/LSD_data_%s_%s_Tmrca.%s"%(data_set, rate_type, fmt))



