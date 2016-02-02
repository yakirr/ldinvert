from __future__ import print_function, division
import numpy as np
import argparse
import pickle
from estimator import Estimator
from primitives import Dataset, GenomicSubset, SnpSubset
# from sparse import ldmatrix
from sparse.blockdiag import BlockDiag


class Truth(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--ld_bandwidth', type=float, required=True,
            help='the maximal ld bandwidth to allow, in Morgans')
    parser.add_argument('--region', type=str, required=True,
            help='the name of the subset of the genome whose heritability should be \
                    analyzed. these files are in data/genome_subsets')

    def readable_name(self):
        return 'Truth,A={},ld_band={}'.format(
                self.params.region,
                self.params.ld_bandwidth)

    def preprocessing_foldername(self):
        return 'pre.truthcovariance.A={}.ldbandwidth={}'.format(
                self.params.region,
                self.params.ld_bandwidth)

    def RA_file(self, mode='rb'):
        return open(self.path_to_preprocessed_data() + 'RA.bd', mode)

    def preprocess(self):
        gs = GenomicSubset(self.params.region)
        ss = SnpSubset(self.refpanel, bedtool=gs.bedtool)
        RA = BlockDiag.ld_matrix(self.refpanel, ss.irs, self.params.ld_bandwidth)
        pickle.dump(RA, self.RA_file(mode='wb'), 2)

    def run(self, beta_num, sim):
        RA = pickle.load(self.RA_file())
        beta = pickle.load(sim.beta_file(beta_num))
        beta = BlockDiag.from_big1darray(beta, RA.irs)
        results = [beta.dot(RA.dot(beta))]
        print(results[-1])
        return results


if __name__ == '__main__':
    import primitives.genome, primitives.simulation, primitives.dataset
    reload(primitives.genome)
    reload(primitives.simulation)
    reload(primitives.dataset)
    from primitives.genome import GenomicSubset, SnpSubset
    from primitives.simulation import SumstatSimulation
    from primitives.dataset import Dataset
    est = Truth(refpanel='tinyGERA', region='tiny', ld_bandwidth=0.1)
    sim = SumstatSimulation('tinyGERA.tiny_inf')
    # est.preprocess()
    est.run(1, sim)
