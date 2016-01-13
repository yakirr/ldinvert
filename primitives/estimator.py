from __future__ import print_function, division
import argparse
import os, inspect
import itertools
from pyutils import fs, bsub
from simulation import SumstatSimulation
import methods
import paths

class Estimator(object):
    def __init__(self, method_name, Nref, **kwargs):
        self.Nref = Nref
        self.method = methods.find_method(method_name)
        if 'command_line_params' in kwargs:
            self.parse_command_line_params(kwargs['command_line_params'])
        else:
            self.params = kwargs

    def parse_command_line_params(self, command_line_params):
        self.params = self.method.parser.parse_known_args(
                command_line_params)[0].__dict__

    def command_line_params(self):
        return list(itertools.chain(*[['--' + n, str(v)] for n, v in self.params.items()]))

    def params_str(self):
        l = [n + '=' + str(v) for n, v in self.params.items()]
        return ','.join(l)

    def __str__(self):
        return self.method.__name__ + ',' + str(self.Nref) + ',' + str(self.params)

    # Preprocessing code
    def path_to_preprocessed_data(self, sim, create=True):
        path = sim.path_to_refpanel() + 'preprocessing.Nref=' + str(self.Nref) + '.' + \
                self.method.required_preprocessing + '/'
        if create:
            fs.makedir(path)
        return path

    def preprocess_job_name(self, sim):
        return 'preprocess-{}-{}-{}'.format(self.method.__name__, self.Nref, sim.name)

    def preprocess_submitted(self, sim):
        return os.path.isfile(self.path_to_preprocessed_data(sim) + '.submitted')

    def declare_preprocess_submitted(self, sim):
        f = open(self.path_to_preprocessed_data(sim) + '.submitted', 'w')
        f.close()

    def submit_preprocess(self, sim):
        if not self.preprocess_submitted(sim):
            print(str(self), 'pre-processing', sim.name)
            my_args = ['--sim_name', sim.name,
                    '--method_name', self.method.__name__,
                    '--Nref', str(self.Nref),
                    'preprocess'] + \
                    self.command_line_params()
            outfilepath = self.path_to_preprocessed_data(sim) + '.preprocessing.out'
            bsub.submit(
                    ['python', '-u', paths.code + 'methods/estimator.py'] + my_args,
                    outfilepath,
                    jobname=self.preprocess_job_name(sim),
                    memory_GB=16)
            self.declare_preprocess_submitted(sim)
        else:
            print(str(self), ': pre-processing unnecessary for', sim.name)

    def do_preprocessing(self, sim):
        self.method.preprocess(sim, self.Nref)

    # Running code
    def path_to_results(self, sim, beta_num):
        return sim.path_to_beta(beta_num, create=False) + 'results.' + \
                'Nref=' + str(self.Nref) + '.' + \
                self.method.__name__ + '.' + self.params_str()

    def run_job_name(self, sim):
        return 'run-{}-{}-{}-{}[1-{}]'.format(
                self.method.__name__, self.params_str(), self.Nref, sim.name, sim.num_betas)

    def submit_runs(self, sim):
        print('\n' + str(self), 'submitting', sim.name)
        my_args = ['--sim_name', sim.name,
                '--method_name', self.method.__name__,
                '--Nref', str(self.Nref),
                'run',
                '--beta_num', '$LSB_JOBINDEX'] + \
                self.command_line_params()
        outfilepath = self.path_to_results(sim, '%I')
        bsub.submit(
                ['python', '-u', paths.code + 'methods/estimator.py'] + my_args,
                outfilepath,
                jobname=self.run_job_name(sim),
                memory_GB=16)

    def run_on_beta(self, beta_num, sim):
        self.method.run(beta_num, sim, self.Nref, **self.params)


def preprocess(est, sim, args):
    est.do_preprocessing(sim)

def run_on_beta(est, sim, args):
    est.run_on_beta(args.beta_num, sim)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sim_name', type=str, required=True,
            help='the name of the simulation for which we are preprocessing')
    parser.add_argument('--method_name', type=str, required=True,
            help='the name of the method for which we are preprocessing')
    parser.add_argument('--Nref', type=int, required=True,
            help='the sample size of the reference panel to allow for the preprocessing')

    subparsers = parser.add_subparsers()
    subparser_preprocess = subparsers.add_parser('preprocess')
    subparser_preprocess.set_defaults(_func=preprocess)
    subparser_run = subparsers.add_parser('run')
    subparser_run.set_defaults(_func=run_on_beta)
    subparser_run.add_argument('--beta_num', type=int, required=True,
            help='the 1-based index of the beta on which the estimator should be run')

    # construct estimator and simulation and do the appropriate action
    args, remaining = parser.parse_known_args()
    est = Estimator(args.method_name, args.Nref, command_line_params=remaining)
    sim = SumstatSimulation(args.sim_name)
    args._func(est, sim, args)
