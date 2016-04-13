from __future__ import print_function, division
import argparse
from estimator import Estimator


class ACorFE(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--annot_chr', type=str, required=True,
            help='path to a set of annot files to which we are computing covariance')
    parser.add_argument('--annot_name', type=str, required=True,
            help='short name to describe the annotation')

    def readable_name(self):
        return 'acorfe,A={}'.format(
                self.params.annot_name)
    def preprocessing_foldername(self):
        return 'acor.preprocess'

    def preprocess(self):
        pass

    def run(self, beta_num, sim):
        print('acor_fe is running on', sim.name, beta_num, 'with')
        print(self.params)

class TruthCorFE(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--annot_chr', type=str, required=True,
            help='path to a set of annot files to which we are computing covariance')
    parser.add_argument('--annot_name', type=str, required=True,
            help='short name to describe the annotation')

    def readable_name(self):
        return 'acorfe,A={}'.format(
                self.params.annot_name)
    def preprocessing_foldername(self):
        return 'acor.preprocess'

    def preprocess(self):
        pass

    def run(self, beta_num, sim):
        print('truthcorfe is running on', sim.name, beta_num, 'with')
        print(self.params)

class TruthCorRE(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--annot_chr', type=str, required=True,
            help='path to a set of annot files to which we are computing covariance')
    parser.add_argument('--annot_name', type=str, required=True,
            help='short name to describe the annotation')

    def readable_name(self):
        return 'acorfe,A={}'.format(
                self.params.annot_name)
    def preprocessing_foldername(self):
        return 'acor.preprocess'

    def preprocess(self):
        pass

    def run(self, beta_num, sim):
        return sim.


if __name__ == '__main__':
    import primitives.genome, primitives.simulation, primitives.dataset
    from primitives.genome import GenomicSubset, SnpSubset
    from primitives.simulation import SumstatSimulation
    from primitives.dataset import Dataset
    est = TruthCorFE(refpanel='UK10Khg19', annot_name='mock',
            annot_chr='/groups/price/yakir/data/signed_annot/1000G_EUR_Phase1/mock')
    sim = SumstatSimulation('test')
    # est.preprocess()
    est.run(1, sim)
