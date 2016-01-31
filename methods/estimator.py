from __future__ import print_function, division
import abc
import numpy as np
import argparse
import pickle
import os
import itertools
from pyutils import fs, bsub
import paths


#TODO: consider moving the code that manages the submission and the files around that
#   to the estimator_manager.py class
class Estimator(object):
    __metaclass__ = abc.ABCMeta
    def __init__(self, **kwargs):
        if 'command_line_params' in kwargs:
            self.parse_command_line_params(kwargs['command_line_params'])
        else:
            self.params = argparse.Namespace(**kwargs)

    def method(self):
        return self.__class__.__name__

    def parse_command_line_params(self, command_line_params):
        self.params = argparse.Namespace(**self.__class__.parser.parse_known_args(
                command_line_params)[0].__dict__)

    def command_line_params(self):
        return list(itertools.chain(*[
            ['--' + n, str(v)] for n, v in self.params.__dict__.items()
            ]))

    def params_str(self):
        l = [n + '=' + str(v) for n, v in self.params.__dict__.items()]
        return ','.join(sorted(l))

    def __str__(self):
        return self.method() + ',' + str(self.params.__dict__)

    @abc.abstractmethod
    def readable_name(self): pass

    # Preprocessing code
    @abc.abstractmethod
    def preprocessing_folder(self): pass

    def path_to_preprocessed_data(self, sim, create=True):
        path = sim.path_to_refpanel() + self.preprocessing_folder() + '/'
        if create:
            fs.makedir(path)
        return path

    def preprocess_job_name(self, sim):
        return 'preprocess-{}-{}-{}'.format(
                self.method(), self.params_str(), sim.name)

    def preprocess_submitted(self, sim):
        return os.path.isfile(self.path_to_preprocessed_data(sim) + '.submitted')

    def declare_preprocess_submitted(self, sim):
        f = open(self.path_to_preprocessed_data(sim) + '.submitted', 'w')
        f.close()

    def submit_preprocess(self, sim):
        if not self.preprocess_submitted(sim):
            print(str(self), 'pre-processing', sim.name)
            my_args = ['--sim_name', sim.name,
                    '--method_name', self.method(),
                    'preprocess'] + \
                    self.command_line_params()
            outfilepath = self.path_to_preprocessed_data(sim) + '.preprocessing.out'
            bsub.submit(
                    ['python', '-u', paths.code + 'methods/estimator_manager.py'] + my_args,
                    outfilepath,
                    jobname=self.preprocess_job_name(sim),
                    memory_GB=16)
            self.declare_preprocess_submitted(sim)
        else:
            print(str(self), ': pre-processing unnecessary for', sim.name)

    @abc.abstractmethod
    def preprocess(self, sim): pass

    # Running code
    def results_name(self):
        return 'results.' + self.method() + '.' + self.params_str()

    def results_path_stem(self, sim, beta_num):
        return sim.path_to_beta(beta_num, create=False) + self.results_name()

    def run_job_name(self, sim):
        return 'run-{}-{}-{}[1-{}]'.format(
                self.method(), self.params_str(), sim.name, sim.num_betas)

    def submit_runs(self, sim):
        print('\n' + str(self), 'submitting', sim.name)
        my_args = ['--sim_name', sim.name,
                '--method_name', self.method(),
                'run',
                '--beta_num', '$LSB_JOBINDEX'] + \
                self.command_line_params()
        outfilepath = self.results_path_stem(sim, '%I') + '.out'
        if all([os.path.exists(self.results_path_stem(sim, beta_num))
            for beta_num in range(1, sim.num_betas+1)]):
            print('submission unnecessary for', str(self))
        else:
            bsub.submit(
                    ['python', '-u', paths.code + 'methods/estimator_manager.py'] + my_args,
                    outfilepath,
                    jobname=self.run_job_name(sim),
                    memory_GB=16)

    @abc.abstractmethod
    def run(self, beta_num, sim): pass

    def run_and_save_results(self, beta_num, sim):
        results = self.run(beta_num, sim)
        print(results)
        np.savetxt(self.results_path_stem(sim, beta_num), results)

    def results(self, beta_num, sim):
        return np.loadtxt(self.results_path_stem(sim, beta_num))
