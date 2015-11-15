from __future__ import print_function, division
import numpy as np
import pickle
import argparse
import subprocess
from time import time
from pysnptools.util import IntRangeSet
import pyutils.fs as fs
import pyutils.configs
import genome.utils as gutils
import sim.utils as simutils
import hyperparams as hp

# returns an array that is genome_size x num_betas
def sample_beta(p_causal_class, p_causal_rest,
                sigma2_class, sigma2_rest, pathway_intrangeset, num_betas=1):
    # returns an array that is num_SNPs x num_samples
    def point_normal(num_samples, p_causal, sigma2, num_SNPs):
        result = np.zeros((num_SNPs, num_samples))
        causal_mask = np.random.multinomial(
                1,
                [p_causal, 1-p_causal],
                size=(num_samples,num_SNPs))[:,:,0].T.astype(bool)
        num_causal = np.sum(causal_mask)
        result[causal_mask] = np.sqrt(sigma2) * np.random.randn(num_causal)
        return result

    snps_in_pathway = len(pathway_intrangeset)
    genome_effects = point_normal(num_betas, p_causal_rest, sigma2_rest,
                                    hp.dataset.M() - snps_in_pathway)
    pathway_effects = point_normal(num_betas, p_causal_class, sigma2_class, snps_in_pathway)

    result = np.empty((hp.dataset.M(), num_betas))
    result[pathway_intrangeset, :] = pathway_effects
    result[hp.dataset.all_snps()-pathway_intrangeset, :] = genome_effects

    return result

def main(args):
    np.random.seed(args.beta_num)
    hp.load()

    # load relevant outside information about pathway and genome as well as genotypes
    pathway_intrangeset = pickle.load(hp.pathway_file()) & hp.dataset.all_snps()
    pathway_size_SNPs = len(pathway_intrangeset)

    # compute the parameters we'll need to sample beta
    if hp.beta_params.pA > 0:
        sigma2A = hp.beta_params.h2gA / (hp.beta_params.pA * pathway_size_SNPs)
    else:
        sigma2A = 0
    if hp.beta_params.pG > 0:
        sigma2G = hp.beta_params.h2gG / (hp.beta_params.pG *
                                        (hp.dataset.M() - pathway_size_SNPs))
    else:
        sigma2G = 0

    # sample the beta
    beta = sample_beta(
        hp.beta_params.pA,
        hp.beta_params.pG,
        sigma2A,
        sigma2G,
        pathway_intrangeset)[:,0]

    # compute noiseless phenotypes slice by slice
    Y = np.zeros(hp.dataset.N)
    t0 = time()
    for s in gutils.slices(hp.dataset):
        # X will be N x M
        print(int(time() - t0), ': getting genotypes from file. SNPs', s)
        X = simutils.get_standardized_genotypes(s)
        print('computing phenotypes. SNPs', s)
        Y += X.dot(beta[s[0]:s[1]])
        del X

    # normalize the Y and the beta to the desired heritability
    normalization = np.std(Y) / np.sqrt(hp.beta_params.h2gA + hp.beta_params.h2gG)
    if normalization == 0: normalization = 1  # just in case we have some 0s...
    Y /= normalization
    beta /= normalization

    # write the betas and the noiseless phenotypes
    pickle.dump(beta, hp.beta_file(args.beta_num, 'wb'), 2)
    pickle.dump(Y, hp.noiseless_Y_file(args.beta_num, 'wb'), 2)

def submit(args):
    hp.load(printall=False)
    my_args = ['main', '--beta_num', '$LSB_JOBINDEX']
    outfilepath = hp.path_to_results_dir() + '.out.beta.%I.noiseless'
    cmd = pyutils.configs.bsub_command(
            ['python', '-u', hp.paths.code + 'sim/sample_betas.py'] + \
                    my_args + hp.to_command_line(),
            outfilepath,
            jobname='samplebetas[1-20]',
            memory_GB=16)
    print(' '.join(cmd))
    print(outfilepath)
    subprocess.call(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    subparser_main = subparsers.add_parser('main')
    subparser_main.add_argument('--beta_num', type=int, required=True,
            help='the 1-based index of the beta to generate')
    subparser_main.set_defaults(_func=main)

    subparser_submit = subparsers.add_parser('submit')
    subparser_submit.set_defaults(_func=submit)

    args, _ = parser.parse_known_args()
    pyutils.configs.print_vars(args)
    print()

    args._func(args)
