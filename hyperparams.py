from __future__ import print_function, division
import pickle
import argparse
from pysnptools.snpreader import Bed
import pyutils.fs as fs
import pyutils.configs

#######################################################
## Global variables that aren't read from command line
#######################################################
dataset = argparse.Namespace()
dataset.name = 'WT1QC_fewnans'
dataset.reference_genome = 'hg18'
dataset.pathway = '99'

paths = argparse.Namespace()
paths.code = '/groups/price/yakir/py/'
paths.genotypes = '/groups/price/hilary/data/' + dataset.name + '/'
paths.reference = '/groups/price/yakir/data/reference/'
paths.pathways_with_flanks = '/groups/price/yakir/data/GO.1Mb_flanks/' + \
        dataset.pathway + '.' + dataset.name + '/'
paths.sumstats = '/groups/price/yakir/data/sumstats/' + \
        dataset.pathway + '.' + dataset.name + '/'

print('==global variables==')
pyutils.configs.print_vars(dataset)
pyutils.configs.print_vars(paths)
print('====================')

#TODO: change name to globals
#######################################################
## Variables read from command line
#######################################################
params_parser = argparse.ArgumentParser()
params_parser.add_argument('--first_chrom', type=int, default=1,
        help='the smallest-numbered chromosome to include in the analysis')
def chromosomes():
    return range(params.first_chrom, 23)

beta_params_parser = argparse.ArgumentParser()
beta_params_parser.add_argument('--h2gA', type=float, default=0.3,
        help='the heritability to include in the pathway, on average')
beta_params_parser.add_argument('--h2gG', type=float, default=0.7,
        help='the heritability to include in the rest of the genome on average')
beta_params_parser.add_argument('--pA', type=float, default=1,
        help='the probability of a SNP in the pathway being causal')
beta_params_parser.add_argument('--pG', type=float, default=1,
        help='the probability of a SNP in the rest of the genome being causal')

sumstats_params_parser = argparse.ArgumentParser()
sumstats_params_parser.add_argument('--N', type=int, default=14526)

#######################################################
## Functions that determine the file structure of the data
#######################################################
def pathway_file(mode='rb'):
    return open(paths.pathways_with_flanks + 'pathway.regions_to_indexsets', mode)

def pathway_with_flanks_file(mode='rb'):
    return open(paths.pathways_with_flanks + 'merged.regions_to_indexsets', mode)

def genotypes_bed_file(chrnum):
    return Bed(paths.genotypes + 'all.' + str(chrnum))

def covariance_around_pathway_file(mode='rb'):
    return open(paths.pathways_with_flanks + 'merged.covariance.bda', mode)

def ldscores_files(mode='rb'):
    return open(paths.pathways_with_flanks + 'pathway.chrnum_to_ldscores', mode), \
            open(paths.pathways_with_flanks + 'notpathway.chrnum_to_ldscores', mode)


def results_dirname():
    return 'fc=' + str(params.first_chrom) + \
        ',pG=' + str(float(beta_params.pG)) + \
        ',pA=' + str(float(beta_params.pA)) + \
        ',h2g=' + str(float(beta_params.h2gG+beta_params.h2gA)) + \
        ',Eh2gA=' + str(float(beta_params.h2gA))

# path_to_results_dir:  99.WT1QC_fewnans/pG=x,pA=x,h2g=x,Eh2gA=x/
def path_to_results_dir():
    return paths.sumstats + results_dirname() + '/'

def results_file(mode='r'):
    return open(path_to_results_dir() + 'results.tsv', mode=mode)

# path_to_beta_dir: 99.WT1QC_fewnans/pG=x,pA=x,h2g=x,Eh2gA=x/beta.0/
# (contains beta file and noiseless Ys for entire dataset
def path_to_beta_dir(beta_num, create=True):
    path = path_to_results_dir() + 'beta.' + str(beta_num) + '/'
    if create:
        fs.makedir(path)
    return path

def beta_file(beta_num, mode='rb'):
    return open(path_to_beta_dir(beta_num) + 'chrnum_to_beta', mode)

def noiseless_Y_file(beta_num, mode='rb'):
    return open(path_to_beta_dir(beta_num) + 'noiseless_Y.1darray', mode)

# path_to_samplesize_dir:   path_to_beta_dir + N=1000/
#       for many samples, contains: noisy phenotypes, individuals in sample, and sumstats
def path_to_samplesize_dir(beta_num, create=True):
    path = path_to_beta_dir(beta_num) + 'N=' + str(sumstats_params.N) + '/'
    if create:
        fs.makedir(path)
    return path

def noisy_Y_file(beta_num, index, mode='rb'):
    return open(path_to_samplesize_dir(beta_num) + str(index) + '.Y.1darray', mode)

def individuals_file(beta_num, index, mode='rb'):
    return open(path_to_samplesize_dir(beta_num) + str(index) + '.indivs_indices.1darray',
            mode)

def sumstats_file(beta_num, index, mode='rb'):
    return open(path_to_samplesize_dir(beta_num) + str(index) + '.chrnum_to_alphahat', mode)

#######################################################
## Functions for loading and printing the command-line parameters
#######################################################
def load(printall=True):
    global params, beta_params, sumstats_params
    params, _ = params_parser.parse_known_args()
    beta_params, _ = beta_params_parser.parse_known_args()
    sumstats_params, _ = sumstats_params_parser.parse_known_args()

    if printall:
        print('==hyperparameters===')
        pyutils.configs.print_vars(params); print()
        pyutils.configs.print_vars(beta_params); print()
        pyutils.configs.print_vars(sumstats_params)
        print('====================')

def to_command_line():
    return pyutils.configs.namespace_to_commandline(params) + \
        pyutils.configs.namespace_to_commandline(beta_params) + \
        pyutils.configs.namespace_to_commandline(sumstats_params)

