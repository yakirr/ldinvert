from __future__ import print_function, division
import numpy as np
import pickle
import argparse
import subprocess
from pysnptools.util import IntRangeSet
import pyutils.fs as fs
import pyutils.configs
import genome.utils as gutils
import sim.utils as simutils
import hyperparams as hp

# returns an array that is genome_size x num_betas
def sample_betas(num_betas, chrnum, p_causal_class, p_causal_rest,
                sigma2_class, sigma2_rest, p_regions_to_indexsets):
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

    SNPs_on_chrom = gutils.chromosome_length_in_SNPs(chrnum, hp.dataset.name)
    p_indexset_wrt_chrom = gutils.indexset_wrt_chromosome(
            chrnum,
            p_regions_to_indexsets,
            hp.dataset.name)

    SNPs_in_pathway = len(p_indexset_wrt_chrom)
    genome_effects = point_normal(num_betas, p_causal_rest, sigma2_rest,
                                    SNPs_on_chrom - SNPs_in_pathway)
    pathway_effects = point_normal(num_betas, p_causal_class, sigma2_class, SNPs_in_pathway)

    result = np.empty((SNPs_on_chrom, num_betas))
    result[p_indexset_wrt_chrom, :] = pathway_effects
    result[IntRangeSet((0, SNPs_on_chrom))-p_indexset_wrt_chrom, :] = genome_effects

    return result

def main(args):
    np.random.seed(0)
    hp.load()

    # load relevant outside information about pathway and genome as well as genotypes
    p_regions_to_indexsets = pickle.load(
            open(hp.paths.pathways_with_flanks + 'pathway.regions_to_indexsets', 'rb'))
    pathway_size_SNPs = gutils.total_size_snps(p_regions_to_indexsets.values())
    genome_size_SNPs = gutils.snps_in_dataset(hp.dataset.name)

    # compute the parameters we'll need to sample beta
    if hp.beta_params.pA > 0:
        sigma2A = hp.beta_params.h2gA / (hp.beta_params.pA * pathway_size_SNPs)
    else:
        sigma2A = 0
    if hp.beta_params.pG > 0:
        sigma2G = hp.beta_params.h2gG / (hp.beta_params.pG *
                                        (genome_size_SNPs - pathway_size_SNPs))
    else:
        sigma2G = 0

    # create the betas chromosome by chromosome
    chrnum_to_betas = {}
    def generate_betas_for_chr(chrnum):
        print('sampling betas. chrom', chrnum)
        chrnum_to_betas[chrnum] = sample_betas(
                args.num_betas,
                chrnum,
                hp.beta_params.pA,
                hp.beta_params.pG,
                sigma2A,
                sigma2G,
                p_regions_to_indexsets)
    map(generate_betas_for_chr, hp.chromosomes())

    # compute noiseless phenotypes chrom by chrom
    Ys = np.zeros((gutils.sample_size(hp.dataset.name), args.num_betas))
    for chrnum in hp.chromosomes():
        # X will be N x M
        print('getting genotypes from file. chrom', chrnum)
        X = simutils.get_standardized_genotypes(chrnum)
        print('computing phenotypes. chrom', chrnum)
        Ys += X.dot(chrnum_to_betas[chrnum])

    # normalize the Ys and the betas to the desired heritability
    normalization = np.std(Ys, axis=0) / np.sqrt(hp.beta_params.h2gA + hp.beta_params.h2gG)
    normalization[normalization == 0] = 1  # just in case we have some 0s...
    Ys /= normalization
    def normalize_betas_for_chr(chrnum):
        chrnum_to_betas[chrnum] /= normalization
    map(normalize_betas_for_chr, hp.chromosomes())

    # write the betas
    def print_beta(beta_num):
        beta_i = {chrnum:betas[:,beta_num] for chrnum, betas in chrnum_to_betas.items()}
        pickle.dump(beta_i, hp.beta_file(beta_num + 1, 'wb'), 2)
    map(print_beta, range(args.num_betas))

    # write the noiseless phenotypes
    def print_Y(beta_num):
        pickle.dump(Ys[:,beta_num], hp.noiseless_Y_file(beta_num + 1, 'wb'), 2)
    map(print_Y, range(args.num_betas))

def submit(args):
    hp.load(printall=False)
    my_args = ['main', '--num_betas', '20']
    outfilepath = hp.paths.sumstats + 'out.' + hp.results_dirname()
    cmd = pyutils.configs.bsub_command(
            ['python', '-u', hp.paths.code + 'sim/sample_betas.py'] + \
                    my_args + hp.to_command_line(),
            outfilepath)
    print(' '.join(cmd))
    print(outfilepath)
    subprocess.call(cmd)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    subparser_main = subparsers.add_parser('main')
    subparser_main.add_argument('--num_betas', type=int, required=True,
            help='the number of random betas to generate')
    subparser_main.set_defaults(_func=main)

    subparser_submit = subparsers.add_parser('submit')
    subparser_submit.set_defaults(_func=submit)

    args, _ = parser.parse_known_args()
    pyutils.configs.print_vars(args)
    print()

    args._func(args)
