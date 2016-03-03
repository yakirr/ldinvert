from __future__ import print_function, division
import numpy as np
import matplotlib
import argparse
import pickle
from estimator import Estimator
from primitives import Dataset, GenomicSubset, SnpSubset
# from sparse import ldmatrix
from sparse.blockdiag import BlockDiag


class Truth(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--ld_bandwidth', type=int, required=True,
            help='the maximal ld bandwidth to allow, in milli-Morgans')
    parser.add_argument('--region', type=str, required=True,
            help='the name of the subset of the genome whose heritability should be \
                    analyzed. these files are in data/genome_subsets')

    def readable_name(self):
        return 'Truth,A={},ldband={}'.format(
                self.params.region,
                self.params.ld_bandwidth)

    def preprocessing_foldername(self):
        return 'pre.truthcovariance.A={}.ldband={}'.format(
                self.params.region,
                self.params.ld_bandwidth)

    def RA_file(self, mode='rb'):
        return open(self.path_to_preprocessed_data() + 'RA.bd', mode)
    def RA_plotfilename(self):
        return self.path_to_preprocessed_data() + 'RA.png'

    def preprocess(self):
        matplotlib.use('Agg')
        gs = GenomicSubset(self.params.region)
        ss = SnpSubset(self.refpanel, bedtool=gs.bedtool)
        RA = BlockDiag.ld_matrix(self.refpanel, ss.irs.ranges(), self.params.ld_bandwidth / 1000.)
        try: # if the plotting has some error we don't want to not save the stuff
            # RA.plot(ss.irs, filename=self.RA_plotfilename())
            pass
        except:
            pass
        pickle.dump(RA, self.RA_file(mode='wb'), 2)

    def run(self, beta_num, sim):
        RA = pickle.load(self.RA_file())
        beta = pickle.load(sim.beta_file(beta_num))
        beta = BlockDiag.from_big1darray(beta, RA.ranges())
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
    est = Truth(refpanel='tinyGERA', region='tiny', ld_bandwidth=100)
    sim = SumstatSimulation('tinyGERA.tiny_inf')
    # est.preprocess()
    est.run(1, sim)
