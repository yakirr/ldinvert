from __future__ import print_function, division
import argparse
from experiment import Experiment


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp-name',
            help='the name of the json experiment file whose results we wish to verify')
    parser.add_argument('-v', action='store_true', default=False,
            help='print the actual filenames of the missing result files')
    parser.add_argument('-s', action='store_true', default=False,
            help='suppress information corresponding to sim/estimator pairs with no problem')
    args = parser.parse_args()

    print('loading experiment', args.exp_name)
    exp = Experiment(args.exp_name)

    print('verifying results')
    def verify(s):
        print(s.name)
        for est in exp.estimators_and_truth:
            if est.missing_results(s):
                print('XXX', str(est))
                if args.v:
                    print('\n\t'.join(est.missing_results(s)))
                print('')
            else:
                if not args.s:
                    print('vvv', str(est))
                    print('')
        print('')
    map(verify, exp.simulations)

