from __future__ import print_function, division
import argparse
import pandas as pd
from estimator import Estimator
import primitives.annotation as pa
import paths


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


if __name__ == '__main__':
    import sim.metadata as sm
    est = TruthRE(refpanel='1000G3.wim5u', coeff='ANNOT', output='total')
    s = sm.Simulation('test')
    print(est.run(s, 1))
