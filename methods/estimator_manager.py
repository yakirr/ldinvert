from __future__ import print_function, division
import argparse
from pyutils import pretty
from primitives import SumstatSimulation
import methods


def preprocess(est, args):
    est.preprocess()

def run_on_batch(est, args):
    sim = SumstatSimulation(args.sim_name)
    pretty.print_namespace(sim); print()
    print('batch=', args.batch_num)
    print(est)
    est.run_and_save_results(args.batch_num, sim)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--method_name', type=str, required=True,
            help='the name of the method for which we are preprocessing')

    subparsers = parser.add_subparsers()
    subparser_preprocess = subparsers.add_parser('preprocess')
    subparser_preprocess.set_defaults(_func=preprocess)
    subparser_run = subparsers.add_parser('run')
    subparser_run.set_defaults(_func=run_on_batch)
    subparser_run.add_argument('--sim_name', type=str, required=True,
            help='the name of the simulation for which we are preprocessing')
    subparser_run.add_argument('--batch_num', type=int, required=True,
            help='the 1-based index of the batch of betas on which the estimator \
                    should be run')

    # construct estimator and simulation and do the appropriate action
    args, remaining = parser.parse_known_args()
    est = methods.find_method(args.method_name)(command_line_params=remaining)
    args._func(est, args)
