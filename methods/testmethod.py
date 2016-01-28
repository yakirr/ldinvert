from __future__ import print_function, division
import argparse
from estimator import Estimator


class TestMethodTypeA(Estimator):
    parser = argparse.ArgumentParser()
    parser.add_argument('--Nref', type=int, required=True,
            help='the size of the reference panel to use')
    parser.add_argument('--myparam', type=float, required=True,
            help='the value of the magical thing!')

    def preprocessing_folder(self):
        return 'preprocessing.Nref={}.testpreprocess'.format(
                self.params.Nref)

    def preprocess(self, sim):
        print('TestMethod is preprocessing', sim.name, 'with Nref=', self.params.Nref)

    def run(self, beta_num, sim):
        print('TestMethodTypeA is running on', sim.name, beta_num, 'with')
        print(self.params)


class TestMethodTypeB(TestMethodTypeA):
    def run(self, beta_num, sim):
        print('TestMethodTypeB is running on', sim.name, beta_num, 'with')
        print(self.params)
