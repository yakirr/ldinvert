from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import pickle
import hyperparams as hp

def plot_beta(beta_num):
    beta = pickle.load(hp.beta_file(beta_num))
    p = pickle.load(hp.pathway_file())
    p_f = pickle.load(hp.pathway_with_flanks_file())
    ind = np.zeros(max(p_f)+1, dtype=int)
    ind[p_f] = np.arange(max(p_f))
    f = np.zeros(len(p_f))
    f[ind[p]] = 1

    plt.plot(beta[p_f])
    plt.plot(f * max(beta[p_f]) / 2)
    plt.show()

def plot_pathway():
    R = hp.covariance_around_pathway_file()
    #TODO

