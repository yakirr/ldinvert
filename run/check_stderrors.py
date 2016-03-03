from __future__ import print_function, division
import numpy as np
import pickle, argparse
import matplotlib.pyplot as plt
from experiment import Experiment


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp_name', type=str, required=True,
            help='the name of the json experiment file whose results we wish to plot')
    args = parser.parse_args()

    print('loading experiment', args.exp_name)
    exp = Experiment(args.exp_name)

    print('analyzing results')
    def analyze_results_for(sim):
        error_lists = {}
        result_lists = {}
        truth_lists = {}
        for e in exp.estimators:
            print(e.params.pretty_name)
            total_coverage = 0
            avg_stderror = 0
            for beta in range(1, sim.num_betas+1):
                try:
                    truths = exp.truth.results(beta, sim)
                    results = e.results(beta, sim)
                    variances = e.variances(beta, sim)
                except:
                    # print('warning: couldnt find results for', exp, e, beta)
                    continue
                variances = np.maximum(0, variances)
                intervals_L = results - 2*np.sqrt(variances)
                intervals_R = results + 2*np.sqrt(variances)
                coverage = np.count_nonzero((truths <= intervals_R) & (truths >= intervals_L))
                coverage /= len(results)
                avg_stderror += np.mean(np.sqrt(variances))
                total_coverage += coverage
            total_coverage /= sim.num_betas
            avg_stderror /= sim.num_betas
            print('\tavg stderror:', avg_stderror)
            print('\t95% CI coverage probability:', total_coverage)

    map(analyze_results_for, exp.simulations)

