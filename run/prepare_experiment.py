from __future__ import print_function, division
import argparse
from primitives import Experiment


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp_name', type=str, required=True,
            help='the name of the json experiment file to do the preprocessing for')
    args = parser.parse_args()

    exp = Experiment(args.exp_name)

    for s in exp.simulations:
        for e in exp.estimators:
            e.submit_preprocess(s)
