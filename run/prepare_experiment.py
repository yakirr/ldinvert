from __future__ import print_function, division
import argparse
from experiment import Experiment


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp_name', type=str, required=True,
            help='the name of the json experiment file to do the preprocessing for')
    args = parser.parse_args()

    print('loading experiment', args.exp_name)
    exp = Experiment(args.exp_name)

    print('submitting preprocessing')
    for e in exp.estimators_and_truth():
        e.submit_preprocess()
