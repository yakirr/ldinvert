from __future__ import print_function, division
import abc
import numpy as np
import math
import argparse
import pickle
import os
import itertools
from pyutils import fs, bsub, memo
import sim.metadata as sm
import paths


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

    @property
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

    @property
    @memo.memoized
    def refpanel(self):
        return sm.Dataset(self.params.refpanel)

    def __str__(self):
        return self.method + ',' + str(self.params.__dict__)

    @abc.abstractmethod
    def fsid(self): pass # an id string that can be used in filenames

    # Preprocessing code
    @abc.abstractmethod
    def required_files(self, s): pass
    @classmethod
    def exists_or_is_submitted(cls, fname):
        return os.path.isfile(fname) or os.path.isfile(fname+'.submitted')
    @classmethod
    def mark_as_submitted(cls, fname):
        f = open(fname+'.submitted', 'w'); f.close()
    def dependencies_satisfied(self, s):
        return all([Estimator.exists_or_is_submitted(f) for f in self.required_files(s)])
    def declare_preprocess_submitted(self, s):
        for f in self.required_files(s):
            Estimator.mark_as_submitted(f)

    @property
    def preprocess_memoryreq_GB(self):
        return 2

    def submit_preprocess(self, s, debug=False):
        if not self.dependencies_satisfied(s):
            print(str(self), 'pre-processing')
            my_args = ['--method-name', self.method,
                    '--sim-name', s.name,
                    'preprocess'] + \
                    self.command_line_params()
            outfilepath = self.refpanel.path + '.' + s.name + '.' + self.fsid() + \
                    '.preprocessing.out'
            bsub.submit(
                ['python', '-u', paths.code + 'sim/methods/estimator_manager.py'] + my_args,
                outfilepath,
                jobname='preprocess-' + self.fsid() + '-' + s.name,
                memory_GB=self.preprocess_memoryreq_GB,
                debug=debug)
            self.declare_preprocess_submitted(s)
        else:
            print(str(self), ': pre-processing unnecessary')

    @abc.abstractmethod
    def preprocess(self, s): pass

    # Running code
    def results_file(self, s, beta_num):
        return s.beta_folder(beta_num, create=False) + 'results.' + self.fsid()

    def submit_runs(self, s, overwrite=False, debug=False):
        def outfile_path(batch_num):
            path = s.root_folder() + 'logs/'
            fs.makedir(path)
            return path + '{}-batch{}.out'.format(self.fsid(), batch_num)
        def run_job_name():
            return 'run-{}-{}[1-{}]'.format(
                self.fsid(), s.name, self.num_batches(s))

        if all(os.path.exists(self.results_file(s, beta_num))
            for beta_num in range(1, s.num_betas+1)) and not overwrite:
            print('submission unnecessary for', str(self))
            return

        print('\n' + str(self), 'submitting', s.name)
        my_args = ['--method-name', self.method,
                '--sim-name', s.name,
                'run',
                '--batch-num', '$LSB_JOBINDEX'] + \
                self.command_line_params()
        outfilepath = outfile_path('%I')
        bsub.submit(
                ['python', '-u', paths.code + 'sim/methods/estimator_manager.py'] + my_args,
                outfilepath,
                jobname=run_job_name(),
                memory_GB=4,
                debug=debug)

    @abc.abstractmethod
    def run(self, s, beta_num): pass

    def num_batches(self, s):
        return min(10, s.num_betas)
    def batch_size(self, s):
        return int(math.ceil(s.num_betas / self.num_batches(s)))
    # batch_num and beta_nums are 1-indexed to comply with LSF job indices
    def betas_in_batch(self, s, batch_num):
        start = (batch_num - 1)*self.batch_size(s) + 1
        end = min(s.num_betas+1, start + self.batch_size(s))
        return range(start, end)
    def run_and_save_results(self, batch_num, s):
        for beta_num in self.betas_in_batch(s, batch_num):
            print('===beta_num', beta_num, '====')
            results = self.run(s, beta_num)
            results.to_csv(self.results_file(s, beta_num),
                    sep='\t', index=False)

    # def results(self, beta_num, sim):
    #     return np.loadtxt(self.results_path_stem(sim, beta_num))
    # def variances(self, beta_num, sim):
    #     return np.loadtxt(self.variances_path_stem(sim, beta_num))
