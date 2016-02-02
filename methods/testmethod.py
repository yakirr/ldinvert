from __future__ import print_function, division
import argparse
from estimator import Estimator


class TestMethodTypeA(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--myparam', type=float, required=True,
            help='the value of the magical thing!')

    def preprocessing_foldername(self):
        return 'preprocessing.testpreprocess'

    def preprocess(self):
        print('TestMethod is preprocessing with refpanel=', self.params.refpanel)

    def run(self, beta_num, sim):
        print('TestMethodTypeA is running on', sim.name, beta_num, 'with')
        print(self.params)


class TestMethodTypeB(TestMethodTypeA):
    def run(self, beta_num, sim):
        print('TestMethodTypeB is running on', sim.name, beta_num, 'with')
        print(self.params)
