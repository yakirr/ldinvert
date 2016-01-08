from __future__ import print_function, division
import numpy as np
import pickle
from pysnptools.snpreader import Bed
from pybedtools import BedTool, Interval
import sparse.blockdiag as bd
import genome.utils as gutils
import hyperparams as hp
from sklearn import cross_validation
from scipy.stats import multivariate_normal
from collections import defaultdict
from pysnptools.util import IntRangeSet

import pdb

hp.load()
merged_intrangeset = pickle.load(hp.pathway_with_flanks_file(custom_name='50'))
pathway_intrangeset = pickle.load(hp.pathway_file(custom_name='50'))
merged_cov = pickle.load(hp.covariance_around_pathway_file(custom_name='50'))

# sample random set of 14k individuals (this is done for ldscore and for covariance)
np.random.seed(0)
sample_size = 14000
iids = np.random.choice(hp.dataset.N, size=sample_size, replace=False)

def zero_outside_pathway(r, R):
    result = R.copy()
    my_intrangeset = IntRangeSet(r)
    intersection_intrangeset = my_intrangeset & pathway_intrangeset
    mask = np.zeros(len(my_intrangeset), dtype=bool)
    for s in intersection_intrangeset.ranges():
        start = my_intrangeset.index(s[0])
        end = start + s[1] - s[0]
        mask[start:end] = True
    result[~mask] = 0
    result.T[~mask] = 0 # done this way for compatibility with 1d arrays
    return result

def band(r, R):
    my_intrangeset = IntRangeSet(r)
    intersection_intrangeset = my_intrangeset & pathway_intrangeset
    left_bandwidth = min(intersection_intrangeset) - r[0]
    right_bandwidth = r[1] - max(intersection_intrangeset)
    width = r[1] - r[0]
    for i, row in enumerate(R):
        if i + right_bandwidth < len(row):
            row[i+right_bandwidth:] = 0
        if i - left_bandwidth >= 0:
            row[:i-left_bandwidth] = 0

def optimal_lambda_for_range(r):
    range_data_on_disk = hp.dataset.genotypes_bedfile()
    range_data = range_data_on_disk[iids, r[0]:r[1]].read()
    print('read in all snps in region. shape =', range_data.val.shape)

    kf = cross_validation.KFold(n=sample_size, n_folds=2, random_state=0, shuffle=True)
    Lambdas = [0.5, 0.3, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005]
    costs = defaultdict(list)
    matrices = defaultdict(lambda : np.zeros(r[1]-r[0],r[1]-r[0]))
    for train_index, test_index in kf: # |train| > |test|
    # for test_index, train_index in kf: # |test| > |train|
        train_data = range_data[train_index,:].read()
        test_data = range_data[test_index,:].read()
        print(train_data.iid_count, test_data.iid_count)
        M = train_data.sid_count
        train_data.standardize(); train_data.standardize();
        test_data.standardize(); test_data.standardize();
        sample_cov = train_data.val.T.dot(train_data.val) / train_data.iid_count
        test_sample_cov = test_data.val.T.dot(test_data.val) / test_data.iid_count

        R = test_sample_cov
        RA_train = zero_outside_pathway(r, sample_cov)
        RA_test = zero_outside_pathway(r, R)

        def cost_on_pathway(Lambda):
            cov = (sample_cov + Lambda * np.eye(M)) / (1 + Lambda)
            covinv_R = np.linalg.solve(cov, R)
            prod = covinv_R.T.dot(RA_train.dot(covinv_R))
            opnorm = np.linalg.norm(prod - RA_test, ord=2)
            costs[Lambda].append(opnorm)
            print(Lambda, np.mean(costs[Lambda]))
        def cost_on_window(Lambda):
            cov = (sample_cov + Lambda * np.eye(M)) / (1 + Lambda)
            R = test_sample_cov
            prod = np.dot(R, np.linalg.solve(cov, R))
            opnorm = np.linalg.norm(prod - R, ord=2)
            matrices[Lambda] = matrices[Lambda] + (prod - R)
            costs[Lambda].append(opnorm)
            # print(Lambda, np.mean(costs[Lambda]))
            print(Lambda, np.linalg.norm(matrices[Lambda]))
        def gaussian_cost_on_window(Lambda):
            cov = (sample_cov + Lambda * np.eye(M)) / (1 + Lambda)
            mean = np.zeros(M)
            dist = multivariate_normal(mean=mean,
                    cov=cov)
            loglikelihood = np.sum(dist.logpdf(test_data.val))
            costs[Lambda].append(loglikelihood)
            print(Lambda, np.mean(costs[Lambda]))

        map(cost_on_window, Lambdas)

it = merged_intrangeset.ranges()
map(optimal_lambda_for_range, it)



# def est_covariance(training_iids):
#     covariance_matrices = {}
#     for r in merged_intrangeset.ranges():
#         print('range', r)

#         chrom_data_on_disk = hp.dataset.genotypes_bedfile()
#         range_data = chrom_data_on_disk[training_iids, r[0]:r[1]].read()
#         print('read in all snps in region. shape =', range_data.val.shape)

#         print('standardizing...')
#         # we need the second call because the first one standardized the variance
#         # before substituting 0s for the nans.
#         range_data.standardize(); range_data.standardize()

#         print('computing covariance...')
#         covariance_matrices[r] = \
#                 range_data.val.T.dot(range_data.val) / range_data.iid_count

#     print('assembling into BlockDiagArray...')
#     return bd.BlockDiagArray(covariance_matrices)

# def assess_lambda(Lambda):
#     def assess_fold():
#         leave_out = np.random.randint(0, 1000)
#         training_iids = iids[range(len(iids)) != leave_out]
#         cov = est_covariance(training_iids)
#         likelihood = cov.
# for Lambda in [0.1, 0.05, 0.0025, 0.00125, 0.000625]:
