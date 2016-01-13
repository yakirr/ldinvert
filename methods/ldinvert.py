from __future__ import print_function, division
import argparse
from primitives import Dataset


class MLE(object):
    required_preprocessing = 'covariance'
    parser = argparse.ArgumentParser()

    @classmethod
    def preprocess(cls, sim, Nref):
        d = Dataset(sim.dataset)
        for s in d.slices(buffer_size=30, slice_size=300):
            print(s)
        #TODO: write the code for computing and storing the banded covariance matrix.
        # there should probably be a class called LdMatrix or something that abstracts away
        # which particular matrix implementation you choose to use.

    @classmethod
    def run(cls, beta_num, sim, Nref, **kwargs):
        print('TestMethodTypeA is running on', sim.name, beta_num, 'with Nref=', Nref, 'and')
        for p, v in kwargs.items():
            print(p, v)
