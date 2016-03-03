from __future__ import print_function, division
import abc
import numpy as np
import math
import argparse
import pickle
import os
import itertools
from pyutils import fs, bsub
from primitives import Dataset
import paths


#TODO: consider moving the code that manages the submission and the files around that
#   to the estimator_manager.py class
class Estimator(object):
    __metaclass__ = abc.ABCMeta

    parser = argparse.ArgumentParser()
    parser.add_argument('--refpanel', type=str, required=True,
            help='the name of the data set to use as reference panel')
    parser.add_argument('--pretty_name', type=str, required=False, default='no_name',
            help='the name of the data set to use as reference panel')

    def __init__(self, **kwargs):
        if 'command_line_params' in kwargs:
            self.parse_command_line_params(kwargs['command_line_params'])
        else:
            self.parse_command_line_params(Estimator.to_command_line(kwargs))
        self.__refpanel = None

    def method(self):
        return self.__class__.__name__

    def parse_command_line_params(self, command_line_params):
        self.params = self.__class__.parser.parse_known_args(command_line_params)[0]

    @classmethod
    def to_command_line(cls, dictionary):
        return list(itertools.chain(*[
            ['--' + n, str(v)] for n, v in dictionary.items()
            ]))

    def command_line_params(self):
        return Estimator.to_command_line(self.params.__dict__)

    def __str__(self):
        return self.method() + ',' + str(self.params.__dict__)

    @property
    def refpanel(self):
        if self.__refpanel is None:
            self.__refpanel = Dataset(self.params.refpanel)
        return self.__refpanel

    @abc.abstractmethod
    def readable_name(self): pass

    # Preprocessing code
    @abc.abstractmethod
    def preprocessing_foldername(self): pass

    def path_to_preprocessed_data(self, create=True):
        path = self.refpanel.path + self.preprocessing_foldername() + '/'
        if create:
            fs.makedir(path)
        return path

    def preprocess_job_name(self):
        return 'preprocess-' + self.readable_name()

    def preprocess_submitted(self):
        return os.path.isfile(self.path_to_preprocessed_data() + '.submitted')

    def declare_preprocess_submitted(self):
        f = open(self.path_to_preprocessed_data() + '.submitted', 'w')
        f.close()

    def preprocess_memoryreq_GB(self):
        return 2

    def submit_preprocess(self):
        if not self.preprocess_submitted():
            print(str(self), 'pre-processing')
            my_args = ['--method_name', self.method(),
                    'preprocess'] + \
                    self.command_line_params()
            outfilepath = self.path_to_preprocessed_data() + '.preprocessing.out'
            bsub.submit(
                    ['python', '-u', paths.code + 'methods/estimator_manager.py'] + my_args,
                    outfilepath,
                    jobname=self.preprocess_job_name(),
                    memory_GB=self.preprocess_memoryreq_GB())
            self.declare_preprocess_submitted()
        else:
            print(str(self), ': pre-processing unnecessary')

    @abc.abstractmethod
    def preprocess(self): pass

    # Running code
    def results_name(self):
        return 'results.' + self.readable_name()
    def variances_name(self):
        return self.results_name() + '_var'

    def results_path_stem(self, sim, beta_num):
        return sim.path_to_beta(beta_num, create=False) + self.results_name()
    def variances_path_stem(self, sim, beta_num):
        return sim.path_to_beta(beta_num, create=False) + self.variances_name()

    def outfile_path(self, sim, batch_num):
        path = sim.path() + 'logs/'
        fs.makedir(path)
        return path + '{}.batch.{}.out'.format(self.readable_name(), batch_num)

    def num_batches(self, sim):
        return min(10, sim.num_betas)
    def batch_size(self, sim):
        return int(math.ceil(sim.num_betas / self.num_batches(sim)))

    # batch_num and beta_nums are 1-indexed to comply with LSF job indices
    def betas_in_batch(self, sim, batch_num):
        start = (batch_num - 1)*self.batch_size(sim) + 1
        end = min(sim.num_betas+1, start + self.batch_size(sim))
        return range(start, end)

    def run_job_name(self, sim):
        return 'run-{}-{}[1-{}]'.format(
            self.readable_name(), sim.name, self.num_batches(sim))

    def submit_runs(self, sim, overwrite=False, debug=False):
        #TODO: have it check whether there are more betas or more samples per beta,
        # and then have it decide whether to submit the betas in parallel or in groups
        # probably its easy to have the actual estimator manager be able to accept groups of
        # betas
        print('\n' + str(self), 'submitting', sim.name)
        my_args = ['--method_name', self.method(),
                'run',
                '--sim_name', sim.name,
                '--batch_num', '$LSB_JOBINDEX'] + \
                self.command_line_params()
        outfilepath = self.outfile_path(sim, '%I')
        if all(os.path.exists(self.results_path_stem(sim, beta_num))
            for beta_num in range(1, sim.num_betas+1)) and not overwrite:
            print('submission unnecessary for', str(self))
        else:
            bsub.submit(
                    ['python', '-u', paths.code + 'methods/estimator_manager.py'] + my_args,
                    outfilepath,
                    jobname=self.run_job_name(sim),
                    memory_GB=4,
                    debug=debug)

    @abc.abstractmethod
    def run(self, beta_num, sim): pass

    def run_and_save_results(self, batch_num, sim):
        for beta_num in self.betas_in_batch(sim, batch_num):
            print('===beta_num', beta_num, '====')
            results = np.array(self.run(beta_num, sim))
            print(results)
            if len(results.shape) > 1: # this means we got variance estimates as well
                np.savetxt(self.results_path_stem(sim, beta_num), results[:,0])
                np.savetxt(self.variances_path_stem(sim, beta_num), results[:,1])
            else:
                np.savetxt(self.results_path_stem(sim, beta_num), results)

    def results(self, beta_num, sim):
        return np.loadtxt(self.results_path_stem(sim, beta_num))
    def variances(self, beta_num, sim):
        return np.loadtxt(self.variances_path_stem(sim, beta_num))
