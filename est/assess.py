from __future__ import print_function, division
import numpy as np
import pickle
import os
from collections import defaultdict
import sparse.blockdiag as bd
import genome.utils as gutils
import hyperparams as hp
import est.estimators as estimators

def process_results(truth, results):
    biases, variances = {}, {}
    for m in estimators.methods.keys():
        biases[m] = np.mean(results[m] - truth)
        variances[m] = np.var(results[m])
    return biases, variances

def analyze_sample(beta_num, index):
    chrnum_to_alphahat = pickle.load(hp.sumstats_file(beta_num, index))
    indiv_indices = pickle.load(hp.individuals_file(beta_num, index))
    Y = pickle.load(hp.noisy_Y_file(beta_num, index))
    N = len(Y)
    return {
        m : estimators.methods[m](chrnum_to_alphahat, indiv_indices, Y, N)
        for m in estimators.methods.keys()
        }

def analyze_beta(beta_num, truth):
    index = 1
    results = defaultdict(list)
    while True:
        try:
            test = hp.noisy_Y_file(beta_num, index); test.close()
        except IOError:
            break
        result = analyze_sample(beta_num, index)
        for m in estimators.methods.keys():
            results[m].append(result[m])
        index += 1
    results = {m:np.array(results[m]) for m in estimators.methods.keys()}
    return process_results(truth, results)

def get_truth(beta_num):
    return estimators.truth(
            pickle.load(hp.beta_file(beta_num)))

def analyze(first_beta=1):
    beta_num = first_beta
    truth_allbetas, biases_allbetas, variances_allbetas = [], [], []
    while os.path.exists(hp.path_to_samplesize_dir(beta_num, create=False)):
        print(beta_num)
        truth = get_truth(beta_num)
        biases, variances = analyze_beta(beta_num, truth)
        print(truth, biases, variances, sep='\n\t')
        biases_allbetas.append(biases)
        variances_allbetas.append(variances)
        truth_allbetas.append(truth)
        beta_num += 1
    return truth_allbetas, biases_allbetas, variances_allbetas

def write_results(truth, biases, variances):
    outfile = hp.results_file(mode='w')
    outfile.write('beta_num\ttruth\t')
    outfile.write('\t'.join([
        m + '_bias\t' + m + '_var\t' + m + '_mse'
        for m in estimators.methods.keys()]) + '\n')
    for b in range(len(truth)):
        outfile.write('{}\t{:.4f}'.format(b+1, truth[b]))
        for m in estimators.methods.keys():
            outfile.write('\t{:.4f}\t{:.4f}\t{:.4f}'.format(
                biases[b][m],
                variances[b][m],
                biases[b][m]**2 + variances[b][m]))
        outfile.write('\n')

if __name__ == '__main__':
    hp.load()
    truth, biases, variances = analyze(first_beta=1)
    write_results(truth, biases, variances)




# print('{0:<15s}{1:>15s}{2:>15s}'.format('', 'estimate', 'fixed version'))
# print('{0:<15s}{1:15.5f}{2:15.5f}'.format('avg bias', np.mean(bias_givenbeta), np.mean(fixed_bias_givenbeta)))
# print('{0:<15s}{1:15.5f}{2:15.5f}'.format('avg var', np.mean(variance_givenbeta), np.mean(fixed_variance_givenbeta)))
# print('{0:<15s}{1:15.5f}{2:15.5f}'.format('MSE', np.mean(MSE), np.mean(fixed_MSE)))

