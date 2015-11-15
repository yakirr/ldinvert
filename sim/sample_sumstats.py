from __future__ import print_function, division
import numpy as np
import pickle
import argparse
import subprocess
from time import time
import pyutils.configs
import genome.utils as gutils
import pyutils.iter as it
import sim.utils as sutils
import hyperparams as hp

parser = argparse.ArgumentParser()

def main(args):
    np.random.seed(args.beta_num + args.sample_num * 10000)
    hp.load()

    # read in noiseless phenotypes
    Y = pickle.load(hp.noiseless_Y_file(args.beta_num))

    # choose individuals and create ensemble of Ys
    # indices = np.random.choice(Y.shape[0], size=(args.num_samples, hp.sumstats_params.N))
    indices = np.random.choice(Y.shape[0], size=(hp.sumstats_params.N,))
    Y = Y[indices]

    # compute how much noise to add
    sigma2e = 1 - (hp.beta_params.h2gA + hp.beta_params.h2gG)
    print('adding noise. sigma2e =', sigma2e)
    Y += np.sqrt(sigma2e) * np.random.randn(*Y.shape)

    alphahat = np.zeros(hp.dataset.M())
    t0 = time()
    def compute_sumstats_for_slice(s):
        # X will be N x M
        print(int(time() - t0), ': getting genotypes from file. SNPs', s)
        X = sutils.get_standardized_genotypes(s)

        print('computing sumstats. SNPs', s)
        alphahat[s[0]:s[1]] = X[indices].T.dot(Y) / hp.sumstats_params.N
        del X
    map(compute_sumstats_for_slice, gutils.slices(hp.dataset))

    # write output
    def write_output():
        pickle.dump(indices,
                hp.individuals_file(args.beta_num, args.sample_num, 'wb'), 2)
        pickle.dump(Y,
                hp.noisy_Y_file(args.beta_num, args.sample_num, 'wb'), 2)
        pickle.dump(alphahat,
                hp.sumstats_file(args.beta_num, args.sample_num, 'wb'), 2)
    write_output()

def submit(args):
    hp.load(printall=False)
    def submit_beta(beta_num):
        my_args = ['main',
                '--beta_num', str(beta_num),
                '--sample_num', '$LSB_JOBINDEX']
        outfilepath = hp.path_to_beta_dir(beta_num) + '.out.sumstats.%I'

        cmd = pyutils.configs.bsub_command(
                ['python', '-u', hp.paths.code + 'sim/sample_sumstats.py'] + \
                        my_args + hp.to_command_line(),
                outfilepath,
                jobname='samplesumstats[1-100]',
                memory_GB=16)
        print(' '.join(cmd))
        print(outfilepath)
        subprocess.call(cmd)
    map(submit_beta, range(1,21))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    subparser_main = subparsers.add_parser('main')
    subparser_main.add_argument('--beta_num', type=int, required=True,
            help='the index specifying which beta to use')
    subparser_main.add_argument('--sample_num', type=int, required=True,
            help='the index specifying which independent sample we\'re on for a given beta')
    subparser_main.set_defaults(_func=main)

    subparser_submit = subparsers.add_parser('submit')
    subparser_submit.set_defaults(_func=submit)

    args, _ = parser.parse_known_args()
    pyutils.configs.print_vars(args)
    print()

    args._func(args)
