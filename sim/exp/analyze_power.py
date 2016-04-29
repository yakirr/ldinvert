from __future__ import print_function
import numpy as np
import pandas as pd
from experiment import Experiment
import argparse

def get_power(exp, sig_level):
    def filename(s):
        return '{}{}.power.tsv'.format(exp.results_folder(), s.name)

    for s in exp.simulations:
        print(s.name)
        power_array = []
        for est in exp.estimators:
            results_df = est.results(s)
            # fix p-values if two-sided test pvals were calculated wrong (legacy issue)
            results_df.loc[results_df['P']>1,'P'] = \
                    2 - results_df.loc[results_df['P']>1,'P']
            reject = results_df['P'] < sig_level
            power_array.append((est.params.pretty_name, np.mean(reject)))
        df = pd.DataFrame(columns=['Estimators', 'rejection_prob'], data=power_array)
        print(df)
        print()
        df.to_csv(filename(s), sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp-name', help='name of the json experiment file')
    parser.add_argument('--sig-level', type=float, help='significance level')
    args = parser.parse_args()

    exp = Experiment(args.exp_name)
    get_power(exp, args.sig_level)
