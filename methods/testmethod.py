from __future__ import print_function, division
import argparse


class TestMethodTypeA(object):
    required_preprocessing = 'covariance'
    parser = argparse.ArgumentParser()
    parser.add_argument('--myparam', type=float, required=True,
            help='the value of the magical thing!')

    @classmethod
    def preprocess(cls, sim, Nref):
        print('TestMethod is preprocessing', sim.name, 'with Nref=', Nref)

    @classmethod
    def run(cls, beta_num, sim, Nref, **kwargs):
        print('TestMethodTypeA is running on', sim.name, beta_num, 'with Nref=', Nref, 'and')
        for p, v in kwargs.items():
            print(p, v)


class TestMethodTypeB(object):
    required_preprocessing = 'covariance'
    parser = argparse.ArgumentParser()
    parser.add_argument('--myparam', type=float, required=True,
            help='the value of the magical thing!')

    @classmethod
    def preprocess(cls, sim, Nref):
        TestMethodTypeA.preprocess(sim, Nref)

    @classmethod
    def run(cls, beta_num, sim, Nref, **kwargs):
        print('TestMethodTypeB is running on', sim.name, beta_num, 'with Nref=', Nref, 'and')
        for p, v in kwargs.items():
            print(p, v)
