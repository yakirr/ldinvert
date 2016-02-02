from __future__ import print_function, division
import subprocess
import argparse
import numpy as np
import pandas as pd
from pyutils import bsub
from estimator import Estimator
from primitives import Dataset, GenomicSubset, SnpSubset
import paths
import ldsc.ldscore.regressions


class LDSC(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--ld_bandwidth', type=float, required=True,
            help='the maximal ld bandwidth to allow, in Morgans')
    parser.add_argument('--region', type=str, required=True,
            help='the name of the subset of the genome whose heritability should be \
                analyzed. these files are in data/genome_subsets')

    def readable_name(self):
        return 'LDSC,A={},ref={},ldband={}'.format(
                self.params.region,
                self.params.refpanel,
                self.params.ld_bandwidth)

    def preprocessing_foldername(self):
        return 'pre.ldsc.A={}.ldbandwidth={}'.format(
                self.params.region,
                self.params.ld_bandwidth)

    def annotation_filename(self, chrnum):
        return self.path_to_preprocessed_data() + \
                self.params.region + '.' + str(chrnum) + '.annot.gz'

    def l2_filestem(self, chrnum):
        return self.path_to_preprocessed_data() + self.params.region + '.' + str(chrnum)

    def ref_ld(self):
        result = []
        for chrnum in self.refpanel.chromosomes():
            df = pd.read_csv(self.l2_filestem(chrnum) + '.l2.ldscore.gz', sep='\t',
                usecols = [self.params.region + 'L2', 'OTHERL2'])
            result.append(np.array(df))
        return np.concatenate(result)

    def preprocess(self):
        gs = GenomicSubset(self.params.region)

        # create the annotation file
        for chrnum in self.refpanel.chromosomes():
            d = Dataset(self.params.refpanel, chrnum=chrnum)
            ss = SnpSubset(d, gs.restricted_to_chrom_bedtool(chrnum))
            SnpSubset.print_subsets(self.annotation_filename(chrnum),
                    [ss], [self.params.region], add_other=True)

        # create the ldscores file
        for chrnum in self.refpanel.chromosomes():
            d = Dataset(self.params.refpanel, chrnum=chrnum)
            ldscores_command = [
                    'python', '-u', paths.foreign + 'ldsc/ldsc.py',
                    '--l2',
                    '--ld-wind-cm', str(self.params.ld_bandwidth),
                    '--bfile', d.genotypes_bedfile().filename,
                    '--annot', self.annotation_filename(chrnum),
                    '--out', self.l2_filestem(chrnum)]
            print(' '.join(ldscores_command))
            outfilepath = self.l2_filestem(chrnum) + '.bsub_out'
            bsub.submit(
                    ldscores_command,
                    outfilepath,
                    jobname=self.preprocessing_foldername()+',chr='+str(chrnum))

    def run(self, beta_num, sim):
        d = Dataset(sim.dataset)
        gs = GenomicSubset(self.params.region)
        ss = SnpSubset(d, bedtool=gs.bedtool)

        ref_ldscores = self.ref_ld()
        # needs to change when baseline model/untypedsnps
        w_ld = np.sum(ref_ldscores, axis=1).reshape((d.M,1))
        N = np.ones((d.M, 1)) * d.N
        M_annot = np.array([[ss.num_snps(), d.M - ss.num_snps()]])

        print(('ref_ldscores shape:{}\nw_ld shape:{}\nN shape:{}\n' + \
                'M_annot shape:{}').format(
                    ref_ldscores.shape,
                    w_ld.shape,
                    N.shape,
                    M_annot.shape))

        results = []
        for alphahat in sim.sumstats_files(beta_num):
            hsqhat = ldsc.ldscore.regressions.Hsq(
                    alphahat.reshape((d.M,1)),
                    ref_ldscores,
                    w_ld,
                    N,
                    M_annot,
                    intercept=1)
            results.append(hsqhat.coef[0] * ss.num_snps())
            print(len(results), results[-1])

        return results


if __name__ == '__main__':
    import primitives.genome, primitives.simulation, primitives.dataset
    reload(primitives.genome)
    reload(primitives.simulation)
    reload(primitives.dataset)
    from primitives.genome import GenomicSubset, SnpSubset
    from primitives.simulation import SumstatSimulation
    from primitives.dataset import Dataset
    est = LDSC(refpanel='tinyGERA_ref2k', region='tiny', ld_bandwidth=0.01)
    sim = SumstatSimulation('tinyGERA.tiny_inf')
    est.preprocess()
    # est.run(1, sim)

