from __future__ import print_function, division
import subprocess, os
import argparse
import numpy as np
import pandas as pd
from pyutils import bsub, fs
from estimator import Estimator
from primitives import Dataset, GenomicSubset, SnpSubset
import paths
import ldsc.ldscore.regressions


class LDSC(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--ld_window', type=int, required=True,
            help='the maximal ld window to consider, in milliMorgans')
    parser.add_argument('--baseline', type=int, required=True,
            help='1 to use the baseline model, 0 to not')
    parser.add_argument('--constrain_intercept', type=int, required=False, default=1,
            help='1 to constrain the intercept, 0 to not')
    parser.add_argument('--region', type=str, required=True,
            help='the name of the subset of the genome whose heritability should be \
                analyzed. these files are in data/genome_subsets')

    baseline_model_regions = ['whole_genome'] + [
        'ldsc_baseline/' + r.strip()
        for r in open(
            paths.genome_subsets + 'ldsc_baseline/baseline_regions.txt').readlines()]

    def readable_name(self):
        return 'LDSC,A={},ref={},ldwindow={},baseline={},constrained={}'.format(
                self.params.region,
                self.params.refpanel,
                self.params.ld_window,
                self.params.baseline,
                self.params.constrain_intercept)

    def baseline_path(self):
        return self.refpanel.path + 'pre.ldsc.baseline.ldwindow={}/'.format(
                self.params.ld_window)

    def baseline_filename(self, chrnum):
        return self.baseline_path() + 'baseline.{}.annot.gz'.format(chrnum)

    def baseline_l2_filestem(self, chrnum):
        return self.baseline_path() + 'baseline.{}'.format(chrnum)

    def baseline_preprocessing_in_progress(self):
        return os.path.exists(self.baseline_path())

    def declare_baseline_preprocessing_in_progress(self):
        fs.makedir(self.baseline_path())

    def preprocessing_foldername(self):
        return 'pre.ldsc.A={}.ldwindow={}'.format(
                self.params.region,
                self.params.ld_window)

    def annotation_filename(self, chrnum):
        return self.path_to_preprocessed_data() + \
                self.params.region + '.' + str(chrnum) + '.annot.gz'

    def annotation_l2_filestem(self, chrnum):
        return self.path_to_preprocessed_data() + self.params.region + '.' + str(chrnum)

    # returns tuple of (ld scores (M x n_annot), weights (M x 1), number of snps (1 x n_annot)
    def ld_score_info(self):
        cols = [self.params.region+'L2'] + (['OTHERL2'] if not self.params.baseline else [])
        ldscores = []
        M_annot = np.zeros((len(cols) +
            (0 if not self.params.baseline else len(LDSC.baseline_model_regions)),))

        # ld scores and number of snps in each category
        for chrnum in self.refpanel.chromosomes():
            df = pd.read_csv(self.annotation_l2_filestem(chrnum) + '.l2.ldscore.gz', sep='\t',
                    usecols=cols)
            M_annot_chr = np.loadtxt(self.annotation_l2_filestem(chrnum)+'.l2.M')[:len(cols)]

            if self.params.baseline:
                df_baseline = pd.read_csv(
                        self.baseline_l2_filestem(chrnum) + '.l2.ldscore.gz', sep='\t',
                        usecols=[r+'L2' for r in LDSC.baseline_model_regions])
                M_annot_baseline_chr = np.loadtxt(
                        self.baseline_l2_filestem(chrnum)+'.l2.M')
                M_annot_chr = np.concatenate([M_annot_chr, M_annot_baseline_chr])
                df = pd.concat([df, df_baseline], axis=1)

            ldscores.append(np.array(df))
            M_annot += M_annot_chr
        ldscores = np.concatenate(ldscores)

        # weights
        if self.params.baseline:
            w_ld = ldscores[:,1].reshape((-1, 1))
        else:
            w_ld = np.sum(ldscores, axis=1).reshape((-1,1))

        return ldscores, w_ld, M_annot.reshape((1,-1))

    def create_baseline_model(self):
        gss = [GenomicSubset(region) for region in LDSC.baseline_model_regions]

        # create the annotation file
        for chrnum in self.refpanel.chromosomes():
            print('creating baseline annot file for chr', chrnum)
            d = Dataset(self.params.refpanel, chrnum=chrnum)
            sss = [SnpSubset(d, gs.restricted_to_chrom_bedtool(chrnum)) for gs in gss]
            SnpSubset.print_subsets(self.baseline_filename(chrnum),
                    sss, LDSC.baseline_model_regions)

        # create the ldscores file
        for chrnum in self.refpanel.chromosomes():
            d = Dataset(self.params.refpanel, chrnum=chrnum)
            ldscores_command = [
                    'python', '-u', paths.foreign + 'ldsc/ldsc.py',
                    '--l2',
                    '--ld-wind-cm', str(self.params.ld_window / 1000.),
                    '--bfile', d.genotypes_bedfile.filename,
                    '--annot', self.baseline_filename(chrnum),
                    '--out', self.baseline_l2_filestem(chrnum)]
            print(' '.join(ldscores_command))
            outfilepath = self.baseline_l2_filestem(chrnum) + '.bsub_out'
            bsub.submit(
                    ldscores_command,
                    outfilepath,
                    jobname='baseline,chr='+str(chrnum))

    def overlap_vector(self):
        counts = np.zeros((
            2 + (len(LDSC.baseline_model_regions)-1 if self.params.baseline else 0),))
        for chrnum in self.refpanel.chromosomes():
            annot = pd.read_csv(self.annotation_filename(chrnum), delim_whitespace=True,
                    compression='gzip')
            annot = annot.ix[:,4:]
            if self.params.baseline:
                annot = annot.ix[:,0]
                annot_baseline = pd.read_csv(self.baseline_filename(chrnum),
                        delim_whitespace=True, compression='gzip')
                annot_baseline = annot_baseline.ix[:,4:]
                annot = pd.concat([annot, annot_baseline], axis=1)
            counts += np.dot(annot.ix[:,0], annot)
        return counts

    def preprocess(self):
        if self.params.baseline and not self.baseline_preprocessing_in_progress():
            print('baseline model not found. creating...')
            self.declare_baseline_preprocessing_in_progress()
            self.create_baseline_model()

        print('submitting ld score jobs for annotation of interest')
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
                    '--ld-wind-cm', str(self.params.ld_window / 1000.),
                    '--bfile', d.genotypes_bedfile.filename,
                    '--annot', self.annotation_filename(chrnum),
                    '--out', self.annotation_l2_filestem(chrnum)]
            print(' '.join(ldscores_command))
            outfilepath = self.annotation_l2_filestem(chrnum) + '.bsub_out'
            bsub.submit(
                    ldscores_command,
                    outfilepath,
                    jobname=self.preprocessing_foldername()+',chr='+str(chrnum))

    def run(self, beta_num, sim):
        print('loading data set and region info')
        d = Dataset(sim.dataset)
        gs = GenomicSubset(self.params.region)
        ss = SnpSubset(d, bedtool=gs.bedtool)

        print('loading ld score info')
        ref_ldscores, w_ld, M_annot = self.ld_score_info()
        N = np.ones((d.M, 1)) * d.N

        print(('ref_ldscores shape:{}\nw_ld shape:{}\nN shape:{}\n' + \
                'M_annot shape:{}').format(
                    ref_ldscores.shape,
                    w_ld.shape,
                    N.shape,
                    M_annot.shape))

        overlaps = self.overlap_vector()
        print('num snps overlapping with each category:', overlaps)
        results = []
        variances = []
        for alphahat in sim.sumstats_files(beta_num):
            alphahat = d.N * alphahat ** 2
            if self.params.constrain_intercept:
                hsqhat = ldsc.ldscore.regressions.Hsq(
                        alphahat.reshape((d.M,1)),
                        ref_ldscores,
                        w_ld,
                        N,
                        M_annot,
                        intercept=1)
            else:
                hsqhat = ldsc.ldscore.regressions.Hsq(
                        alphahat.reshape((d.M,1)),
                        ref_ldscores,
                        w_ld,
                        N,
                        M_annot)
            results.append(hsqhat.coef.dot(overlaps))
            variances.append(overlaps.dot(hsqhat.coef_cov).dot(overlaps))
            print('intercept:', hsqhat.intercept)
            print(len(results), results[-1], variances[-1])

        return np.concatenate([np.array([results]).T, np.array([variances]).T], axis=1)


if __name__ == '__main__':
    import primitives.genome, primitives.simulation, primitives.dataset
    reload(primitives.genome)
    reload(primitives.simulation)
    reload(primitives.dataset)
    from primitives.genome import GenomicSubset, SnpSubset
    from primitives.simulation import SumstatSimulation
    from primitives.dataset import Dataset
    est = LDSC(refpanel='tinyGERA_ref2k', region='tiny', ld_window=10)
    sim = SumstatSimulation('tinyGERA.tiny_inf')
    # est.preprocess()
    # est.run(1, sim)

