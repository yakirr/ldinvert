from __future__ import print_function, division
import pickle
import argparse
from pybedtools import BedTool
from pysnptools.snpreader import Bed
from pysnptools.util import IntRangeSet
import pyutils.fs as fs
import pyutils.configs

#######################################################
## data set definitions
#######################################################
WT1QC_fewnans = argparse.Namespace()
WT1QC_fewnans.name = 'WT1QC_fewnans'
WT1QC_fewnans.reference_genome = 'hg18'
WT1QC_fewnans.path = '/groups/price/hilary/data/WT1QC_fewnans/'
WT1QC_fewnans.genotypes_bedfile = \
        lambda chrnum : Bed(WT1QC_fewnans.path + 'all.' + str(chrnum))
# we use lambda below so that this object isn't loaded when hyperparams is
WT1QC_fewnans.snps_bedtool = lambda : BedTool(WT1QC_fewnans.path + 'all.all_chr.ucscbed')

GERA = argparse.Namespace()
GERA.name = 'GERA'
GERA.reference_genome = 'hg18'
GERA.path = '/groups/price/yakir/data/genotypes/GERA/'
# we use lambda below so that this object isn't loaded when hyperparams is
GERA.genotypes_bedfile = lambda : Bed(GERA.path + 'eur-filtered.all_chr.maf0.01')
GERA.snps_bedtool = lambda: \
    BedTool('/groups/price/yakir/data/genotypes/GERA/eur-filtered.all_chr.maf0.01.ucscbed')

#######################################################
## Global variables that aren't read from command line
#######################################################
dataset = GERA
dataset.N = dataset.genotypes_bedfile().iid_count
dataset.M = lambda : params.trunc_genome if params.trunc_genome > 0 \
                        else dataset.genotypes_bedfile().sid_count
dataset.slice_size = 10000
dataset.all_snps = lambda : IntRangeSet((0, dataset.M()))

pathway = argparse.Namespace()
pathway.name = '99'
pathway.window_size_Mb = 1
pathway.group = str(pathway.window_size_Mb) + 'Mb_flanks'

paths = argparse.Namespace()
paths.code = '/groups/price/yakir/py/'
paths.data = '/groups/price/yakir/data/'
paths.aggregate_results = '/groups/price/yakir/results/'
paths.reference = paths.data + 'reference/'
paths.pathway_details_for = lambda custom_name : paths.data + pathway.group + '/' + \
        custom_name + '.' + dataset.name + '/'
paths.pathway_details = paths.pathway_details_for(pathway.name)
paths.sumstats = '/groups/price/yakir/data/sumstats/' + \
        pathway.name + '.' + dataset.name + '/'
pyutils.fs.makedir(paths.sumstats)

print('==global variables==')
pyutils.configs.print_vars(dataset); print()
pyutils.configs.print_vars(pathway); print()
pyutils.configs.print_vars(paths); print()
print('====================')

#######################################################
## Variables read from command line
#######################################################
params_parser = argparse.ArgumentParser()
params_parser.add_argument('--trunc_genome', type=int, default=-1,
        help='whether to truncate the genome to its first few bp. -1 means full genome.')

beta_params_parser = argparse.ArgumentParser()
beta_params_parser.add_argument('--pG', type=float, default=0.001,
        help='the probability of a SNP in the rest of the genom5 being causal')
beta_params_parser.add_argument('--pA', type=float, default=0.01,
        help='the probability of a SNP in the pathway being causal')
beta_params_parser.add_argument('--h2gG', type=float, default=0.45,
        help='the heritability to include in the rest of the genome on average')
beta_params_parser.add_argument('--h2gA', type=float, default=0.05,
        help='the heritability to include in the pathway, on average')

sumstats_params_parser = argparse.ArgumentParser()
sumstats_params_parser.add_argument('--N', type=int, default=dataset.N)


#######################################################
## Functions that determine the file structure of the data
#######################################################
def pathway_ucscbedfilename():
    return paths.pathway_details + 'pathway.ucscbed'
def pathway_with_flanks_ucscbedfilename():
    return paths.pathway_details + 'merged.ucscbed'
def pathway_flanks_ucscbedfilename():
    return paths.pathway_details + 'flanks.ucscbed'

def pathway_file(custom_name=pathway.name, mode='rb'):
    return open(paths.pathway_details_for(custom_name) + 'pathway.intrangeset', mode)
def pathway_with_flanks_file(custom_name=pathway.name, mode='rb'):
    return open(paths.pathway_details_for(custom_name) + 'merged.intrangeset', mode)
def pathway_flanks_file(custom_name=pathway.name, mode='rb'):
    return open(paths.pathway_details_for(custom_name) + 'flanks.intrangeset', mode)

def covariance_around_pathway_file(custom_name=pathway.name, mode='rb'):
    return open(paths.pathway_details_for(custom_name) + 'merged.covariance.bda', mode)

def ldscores_file_pathway(custom_name=pathway.name, mode='rb'):
    return open(paths.pathway_details_for(custom_name) + 'pathway.ldscores', mode)
def ldscores_file_notpathway(custom_name=pathway.name, mode='rb'):
    return open(paths.pathway_details_for(custom_name) + 'notpathway.ldscores', mode)


def results_dirname():
    return 'pG=' + str(float(beta_params.pG)) + \
        ',pA=' + str(float(beta_params.pA)) + \
        ',h2g=' + str(float(beta_params.h2gG+beta_params.h2gA)) + \
        ',Eh2gA=' + str(float(beta_params.h2gA)) + \
        (',small' if params.trunc_genome > 0 else '')

# path_to_results_dir:  99.WT1QC_fewnans/pG=x,pA=x,h2g=x,Eh2gA=x/
def path_to_results_dir():
    return paths.sumstats + results_dirname() + '/'

def results_file(mode='r'):
    return open(path_to_results_dir() + 'results.tsv', mode=mode)

def results_file_processed(mode='r'):
    return open(path_to_results_dir() + 'results_processed.tsv', mode=mode)

# path_to_beta_dir: 99.WT1QC_fewnans/pG=x,pA=x,h2g=x,Eh2gA=x/beta.0/
# (contains beta file and noiseless Ys for entire dataset
def path_to_beta_dir(beta_num, create=True):
    path = path_to_results_dir() + 'beta.' + str(beta_num) + '/'
    if create:
        fs.makedir(path)
    return path

def beta_file(beta_num, mode='rb'):
    return open(path_to_beta_dir(beta_num) + 'beta', mode)

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
    return open(path_to_samplesize_dir(beta_num) + str(index) + '.alphahat', mode)

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

