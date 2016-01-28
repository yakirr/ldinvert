from __future__ import print_function, division
import argparse
from experiment import Experiment


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp_name', type=str, required=True,
            help='the name of the json experiment file to run')
    args = parser.parse_args()

    exp = Experiment(args.exp_name)

    for s in exp.simulations:
        for e in exp.estimators:
            e.submit_runs(s)
