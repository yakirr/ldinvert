from __future__ import print_function, division
import argparse
import subprocess
import pandas as pd
import scipy.stats as stats
import numpy as np
from pyutils import bsub, pretty
import metadata as sm

# create beta and noiseless phenotype for each chromosome
def create_beta_and_profiles(s, beta_num):
    for chrnum in s.chromosomes:
        # create beta, normalize by maf since genotypes will be centered but not standardized
        print('creating beta for chr', chrnum)
        beta = s.architecture.draw_beta(chrnum)
        maf = s.dataset.frq_df(chrnum)['MAF'].values
        beta['BETA'] /= np.sqrt(maf * (1-maf))
        beta.to_csv(s.beta_file(beta_num, chrnum, mode='w'), index=False, sep='\t')

        # write sparse beta as well (plink works faster if beta is explicitly sparse)
        sparsebeta = beta.loc[beta['BETA'] != 0]
        sparsebeta.to_csv(s.sparsebeta_file(beta_num,chrnum,mode='w'),
                index=False, sep='\t')
        print('norm of beta{}*sqrt(maf(1-maf)) equals {}'.format(
            chrnum, np.linalg.norm(beta['BETA'] * np.sqrt(maf * (1-maf)))**2))
        print('beta{} is length {} and has support of size {}'.format(
            chrnum, len(beta), len(sparsebeta)))

        # call plink to create profile file
        cmd = ['plink',
                '--bfile', s.dataset.bfile(chrnum),
                '--score', s.sparsebeta_filename(beta_num, chrnum), '1', '2', '4',
                'sum',
                'center',
                '--out', s.chr_filestem(beta_num, chrnum)]

        print('executing', ' '.join(cmd))
        subprocess.call(cmd)

# sum the noiseless phenotypes across chromosomes
def make_noiseless_pheno(s, beta_num):
    def get_profile(fname):
        df = pd.read_csv(fname, header=0, delim_whitespace=True)
        df.drop(['PHENO', 'CNT', 'CNT2'], axis=1, inplace=True)
        return df

    print('merging phenotypes. chr', s.chromosomes[0])
    phenotype = get_profile(s.noiselessYchr_filename(beta_num, s.chromosomes[0]))
    phenotype.rename(columns={'SCORESUM' : 'PHENO'}, inplace=True)
    for chrnum in s.chromosomes[1:]:
        print('merging phenotypes. chr', chrnum)
        profile = get_profile(s.noiselessYchr_filename(beta_num, chrnum))
        phenotype = pd.merge(phenotype, profile, how='inner', on=['FID', 'IID'])
        phenotype['PHENO'] += phenotype['SCORESUM']
        phenotype.drop(['SCORESUM'], axis=1, inplace=True)

    phenotype.to_csv(s.noiselessY_filename(beta_num),
            index=False,
            sep='\t')
    return phenotype

# add noise to resulting phenotype
def add_noise_and_save(s, beta_num, phenotype):
    sigma2e = 1-s.h2g
    print('adding noise. sigma2e =', sigma2e)
    phenotype['PHENO'] += np.sqrt(sigma2e) * np.random.randn(len(phenotype))
    phenotype.to_csv(s.noisyY_filename(beta_num),
            index=False,
            sep='\t')

# call plink to compute sumstats
def make_qassoc(s, beta_num):
    # compute one set of sumstats per chromosome
    for chrnum in s.chromosomes:
        print('computing sumstats for chr', chrnum)
        cmd = ['plink',
                '--bfile', s.dataset.bfile(chrnum),
                '--pheno', s.noisyY_filename(beta_num),
                '--allow-no-sex',
                '--assoc',
                '--out', s.chr_filestem(beta_num, chrnum)]
        print('executing', ' '.join(cmd))
        subprocess.call(cmd)

def make_sumstats(s, beta_num):
    # create one large df from the qassoc files.
    print('merging chromosomes of sumstats')
    sumstats = pd.DataFrame()
    for chrnum in s.chromosomes:
        qassoc = pd.read_csv(s.sumstatschr_filename(beta_num, chrnum),
            header=0,
            delim_whitespace=True)
        qassoc['A1'] = s.dataset.bim_df(chrnum)['A1']
        qassoc['A2'] = s.dataset.bim_df(chrnum)['A2']
        sumstats = sumstats.append(qassoc[['SNP', 'A1', 'A2', 'BETA', 'P', 'NMISS']])

    # turn the association statistics into Z scores and match sumstats format
    print('adding sign')
    sumstats['SIGN'] = np.sign(sumstats['BETA'])
    print('computing chisq statistics')
    sumstats['CHISQ'] = stats.chi2.isf(sumstats['P'], 1)
    print('computing z scores')
    sumstats['Z'] = np.sqrt(sumstats['CHISQ']) * sumstats['SIGN']
    sumstats['N'] = sumstats['NMISS']
    sumstats.drop(['BETA', 'P', 'SIGN', 'CHISQ', 'NMISS'], axis=1, inplace=True)

    # save the big sumstats file
    sumstats.to_csv(s.sumstats_file(beta_num, mode='w'), sep='\t', index=False)

def main(args):
    np.random.seed(args.beta_num)
    s = sm.Simulation(args.sim_name)
    create_beta_and_profiles(s, args.beta_num)
    phenotype = make_noiseless_pheno(s, args.beta_num)
    add_noise_and_save(s, args.beta_num, phenotype)
    make_qassoc(s, args.beta_num)
    make_sumstats(s, args.beta_num)

def submit(args):
    s = sm.Simulation(args.sim_name)
    my_args = ['--sim-name', args.sim_name,
            'main',
            '--beta-num', '$LSB_JOBINDEX']
    outfilename = s.beta_folder('%I', create=False) + '.sim_sumstats.out'
    bsub.submit(['python', '-u', __file__] + my_args,
            outfilename,
            jobname='simsumstats[1-{}]'.format(s.num_betas),
            debug=args.debug)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sim-name', type=str, required=True,
            help='name of the set of the json file of simulation parameters, without ext')

    main_parser, submit_parser = bsub.add_main_and_submit(parser, main, submit)

    main_parser.add_argument('--beta-num', type=int, required=True,
            help='index of the beta to simulate. 1-indexed!')
    submit_parser.add_argument('-debug', action='store_true', default=False)

    bsub.choose_parser_and_run(parser)


