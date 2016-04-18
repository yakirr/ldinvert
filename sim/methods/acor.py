from __future__ import print_function, division
import subprocess, os
import argparse
import numpy as np
import pandas as pd
from pyutils import bsub
from estimator import Estimator
import primitives.annotation as pa
import paths


class Acor(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--annot_chr', type=str, required=True,
            help='path to a comma-delimited set of annotgz files, ' + \
                    'not including chromosome number')
    parser.add_argument('--coeff', type=str, required=True,
            help='the name of the coefficient to report')
    parser.add_argument('--kind', type=str, required=True,
            help='re for random effects, fe for fixed')

    def fsid(self):
        return 'Acor.{},coeff={},A={}'.format(
                self.params.kind,
                self.params.coeff,
                self.params.annot_chr.replace('/','_'))

    # NOTE: the function below assumes that the ldscores for the reference panel have
    # already been computed
    def required_files(self, s):
        annots = [pa.Annotation(paths.annotations + aname)
                for aname in self.params.annot_chr.split(',')]
        return [a.conv_filename(c) for c in s.chromosomes for a in annots]

    def preprocess(self, s):
        print('Acor is preprocessing', s.name,
                'with refpanel=', self.params.refpanel)
        print(self.params)

        annots = [pa.Annotation(paths.annotations + aname)
                for aname in self.params.annot_chr.split(',')]

        for a in annots:
            print('preprocessing', a.filestem())
            cmd = [
                    'python', '-u', paths.code + 'acor/acor.py',
                    '--annot-chr', a.stem_chr,
                    '--bfile-chr', self.refpanel.bfile_chr,
                    'conv',
                    '--chroms', ' '.join(str(c) for c in s.chromosomes)]
            print(' '.join(cmd))
            subprocess.call(cmd)

    def run(self, s, beta_num):
        print('Acor is running', s.name, 'on beta', beta_num,
                'with refpanel=', self.params.refpanel)
        print(self.params)
        annots = [pa.Annotation(paths.annotations + aname)
                for aname in self.params.annot_chr.split(',')]

        cmd = [
                'python', '-u', paths.code + 'acor/acor.py',
                '--annot-chr', ' '.join([a.stem_chr for a in annots]),
                '--bfile-chr', self.refpanel.bfile_chr,
                'cor',
                '--ldscores-chr', self.refpanel.bfile_chr,
                '--sumstats', s.sumstats_filename(beta_num),
                '--out', self.results_file(s, beta_num),
                self.params.kind,
                '--chroms', ' '.join(str(c) for c in s.chromosomes)]
        print(' '.join(cmd))
        subprocess.call(cmd)

        acorresults = pd.read_csv(self.results_file(s, beta_num)+'.results',
                delim_whitespace=True,
                header=0)
        rowindex = np.where(np.concatenate([
            list(a.names(s.chromosomes[-1])) for a in annots]) == self.params.coeff)[0]
        estimate = acorresults['MU_EST'][rowindex]
        stderr = acorresults['MU_STDERR'][rowindex]
        pval = acorresults['MU_P'][rowindex]

        return pd.DataFrame(columns=['ESTIMATE', 'STDERR', 'P'],
                data=[[estimate, stderr, pval]])


if __name__ == '__main__':
    import sim.metadata as sm
    est = Acor(refpanel='1000G3.wim5u',
            annot_chr='1000G3.wim5u/mock/',
            kind='re',
            coeff='mock')
    s = sm.Simulation('test')
    # est.preprocess(s)
    est.run(s, 1)

