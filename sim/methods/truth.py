from __future__ import print_function, division
import argparse
import pandas as pd
from estimator import Estimator
import primitives.annotation as pa
import paths
import numpy as np
import pandas as pd

class TruthRE(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--coeff', type=str, required=True,
            help='the coefficient to report')
    parser.add_argument('--output', type=str, required=True,
            help='mean: output mu; meansq: output mu^2; total: output mu^2+sigma^2')

    def fsid(self):
        return 'truthRE-{}'.format(
                self.params.coeff)

    def required_files(self, s):
        return []

    def preprocess(self, s):
        print('nothing to pre-process')

    def run(self, s, beta_num):
        print('TruthRE is running', s.name, 'on beta', beta_num,
                'with refpanel=', self.params.refpanel)
        print(self.params)

        #TODO: below we currently assume only one variance effect. Maybe that's okay though
        mean_effects, var_effects = s.architecture.params()
        if self.params.output == 'mean':
            result = mean_effects[self.params.coeff]
        elif self.params.output == 'meansq':
            result = mean_effects[self.params.coeff]**2
        else: # output total per-snp heritability
            result = mean_effects[self.params.coeff]**2 + \
                    var_effects['ALL']
        return pd.DataFrame(columns=['ESTIMATE', 'STDERR', 'P'],
                data=[[result, 0, 0]])

class TruthFE(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--coeff', type=str, required=True,
            help='the coefficient to report')
    parser.add_argument('--annot_chr', type=str, required=True,
            help='path to annotgz files, not incl chromosome number and extension')

    def fsid(self):
        return 'truthFE-{}'.format(
                self.params.coeff)

    def required_files(self, s):
        return []

    def preprocess(self, s):
        print('nothing to pre-process')

    def run(self, s, beta_num):
        print('TruthFRE is running', s.name, 'on beta', beta_num,
                'with refpanel=', self.params.refpanel)
        print(self.params)
        a = pa.Annotation(paths.annotations + self.params.annot_chr)

        beta = pd.concat(
            [pd.read_csv(s.beta_file(beta_num, chrom), sep = "\t")
                for chrom in s.chromosomes],
            axis = 0)['BETA']
        Rv = pd.concat(
            [pd.read_csv(a.conv_filename(chrom), sep = "\t")
                for chrom in s.chromosomes],
            axis = 0)[self.params.coeff + '.conv1']
        v = pd.concat(
            [a.sannot_df(chrom)
                for chrom in s.chromosomes],
            axis = 0)[self.params.coeff]
        result = np.dot(beta, Rv) / np.dot(v, Rv)
        return pd.DataFrame(columns=['ESTIMATE', 'STDERR', 'P'],
                data=[[result, 0, 0]])



if __name__ == '__main__':
    import sim.metadata as sm
    # est = TruthRE(refpanel='1000G3.wim5u', coeff='ANNOT', output='total')
    est = TruthFE(refpanel='1000G3.wim5u', coeff='ANNOT', annot_chr='1000G3.wim5u/mock/')
    s = sm.Simulation('smalltest')
    # print(est.run(s, 1))

    print(est.results(s))

