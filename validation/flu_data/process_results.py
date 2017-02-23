import pandas
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
plt.ion()

def read_lsd_dataset(fname):

    lsd_cols = ['File', 'N', 'Tmrca_sim', 'mu_sim', 'Runtime', 'objective']
    lsd_df = pandas.read_csv(fname, names=lsd_cols)
    return lsd_df

def read_treetime_dataset(fname):

    cols = ['File', 'N', "Tmrca_sim", "mu_sim", "R2_leaves", "R2_internal", "Runtime"]
    df = pandas.read_csv(fname, names=cols)

    return df



if __name__ == "__main__":

    tt_df = read_treetime_dataset('./H3N2_HA_1980_2015_RES_treetime.csv')
    lsd_df = read_lsd_dataset('./H3N2_HA_1980_2015_RES_lsd.csv')

    Ns = np.unique(np.concatenate((tt_df["N"], lsd_df['N'])))

    Tmrca_mean_tt = []
    Tmrca_err_tt = []
    Mu_mean_tt = []
    Mu_err_tt = []

    Tmrca_mean_lsd = []
    Tmrca_err_lsd = []
    Mu_mean_lsd = []
    Mu_err_lsd = []

    Runtime_mean_lsd = []
    Runtime_mean_tt = []
    Runtime_err_lsd = []
    Runtime_err_tt = []


    for N in Ns:

        Nidx_tt = tt_df["N"] == N
        Tmrca_mean_tt.append((N, tt_df[Nidx_tt]["Tmrca_sim"].mean()))
        Tmrca_err_tt.append((N, tt_df[Nidx_tt]["Tmrca_sim"].std()))
        Mu_mean_tt.append((N, tt_df[Nidx_tt]["mu_sim"].mean()))
        Mu_err_tt.append((N, tt_df[Nidx_tt]["mu_sim"].std()))
        Runtime_mean_tt.append((N, tt_df[Nidx_tt]["Runtime"].mean()))
        Runtime_err_tt .append((N, tt_df[Nidx_tt]["Runtime"].std()))

        Nidx_lsd = lsd_df["N"] == N
        Tmrca_mean_lsd.append((N, lsd_df[Nidx_lsd]["Tmrca_sim"].mean()))
        Tmrca_err_lsd.append((N, lsd_df[Nidx_lsd]["Tmrca_sim"].std()))
        Mu_mean_lsd.append((N, lsd_df[Nidx_lsd]["mu_sim"].mean()))
        Mu_err_lsd.append((N, lsd_df[Nidx_lsd]["mu_sim"].std()))
        Runtime_mean_lsd.append((N, lsd_df[Nidx_lsd]["Runtime"].mean()))
        Runtime_err_lsd .append((N, lsd_df[Nidx_lsd]["Runtime"].std()))


    Tmrca_mean_tt = np.array(Tmrca_mean_tt)
    Tmrca_err_tt = np.array(Tmrca_err_tt)
    Mu_mean_tt = np.array(Mu_mean_tt)
    Mu_err_tt = np.array (Mu_err_tt)

    Tmrca_mean_lsd = np.array(Tmrca_mean_lsd)
    Tmrca_err_lsd = np.array(Tmrca_err_lsd)
    Mu_mean_lsd = np.array(Mu_mean_lsd)
    Mu_err_lsd =  np.array (Mu_err_lsd)

    Runtime_mean_lsd = np.array(Runtime_mean_lsd)
    Runtime_mean_tt  = np.array(Runtime_mean_tt)
    Runtime_err_lsd  = np.array(Runtime_err_lsd)
    Runtime_err_tt   = np.array(Runtime_err_tt)


    plt.figure()
    plt.errorbar(Tmrca_err_tt[:, 0], Tmrca_mean_tt[:, 1], Tmrca_err_tt[:, 1], fmt='o-', label='TreeTime')
    plt.errorbar(Tmrca_err_lsd[:, 0], Tmrca_mean_lsd[:, 1], Tmrca_err_lsd[:, 1], fmt='o-', label='LSD')
    plt.grid()
    plt.legend()
    plt.ylabel("Tmrca")
    plt.title("Estimated Tmrca as function of sample size")
    plt.xlabel("Sample size, MAX # sequences per year")


    plt.figure()
    plt.errorbar(Mu_err_tt[:, 0],  Mu_mean_tt[:, 1],  Mu_err_tt[:, 1], fmt='o-', label='TreeTime')
    plt.errorbar(Mu_err_lsd[:, 0], Mu_mean_lsd[:, 1], Mu_err_lsd[:, 1], fmt='o-', label='LSD')
    plt.grid()
    plt.legend()
    plt.ylabel("Mutation rate")
    plt.title("Estimated Mutation rate as function of sample size")
    plt.xlabel("Sample size, MAX # sequences per year")


    plt.figure()
    plt.errorbar(Runtime_err_tt[:, 0],  Runtime_mean_tt[:, 1],  Runtime_err_tt[:, 1], fmt='o-', label='TreeTime')
    plt.errorbar(Runtime_err_lsd[:, 0], Runtime_mean_lsd[:, 1], Runtime_err_lsd[:, 1], fmt='o-', label='LSD')
    plt.grid()
    plt.legend()
    plt.ylabel("Runtime, arb.u.")
    plt.title("Estimated Runtime as function of sample size")
    plt.xlabel("Sample size, MAX # sequences per year")



