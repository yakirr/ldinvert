from __future__ import print_function, division
import numpy as np
import pickle, argparse, subprocess
from time import time
from pyutils import bsub, pretty
from primitives import SumstatSimulation, Dataset
import paths


def main(args):
    np.random.seed(args.beta_num + args.sample_num * 10000)
    sim = SumstatSimulation(args.sim_name)
    pretty.print_namespace(sim); print()

    # read in noiseless phenotypes
    Y = pickle.load(sim.noiseless_Y_file(args.beta_num))

    # choose individuals and create ensemble of Ys
    indices = np.random.choice(Y.shape[0], size=(sim.sample_size,))
    Y = Y[indices]

    # compute how much noise to add
    sigma2e = 1 - sim.h2g
    print('adding noise. sigma2e =', sigma2e)
    Y += np.sqrt(sigma2e) * np.random.randn(*Y.shape)

    alphahat = np.array([])
    for chrnum in sim.chromosomes:
        print('chr', chrnum)
        d = Dataset(sim.dataset+'.'+str(chrnum))

        # TODO: this should be done only once, not once per chromosome
        # if sim.condition_on_covariates:
        #     print('projecting covariates out of Y')
        #     Y = d.project_out_covariates(Y, covariates=d.covariates[indices])

        t0 = time()
        myalphahat = np.zeros(d.M)
        def compute_sumstats_for_slice(s):
            # X will be N x M
            print(int(time() - t0), ': getting genotypes from file. SNPs', s)
            X = d.get_standardized_genotypes(s)[indices]

            # TODO: projecting out covariates not supported, see above
            # if sim.condition_on_covariates:
            #     print(int(time() - t0), ': projecting out covariates')
            #     X = d.project_out_covariates(X, covariates=d.covariates[indices])

            print(int(time() - t0), ': computing sumstats. SNPs', s)
            myalphahat[s[0]:s[1]] = X.T.dot(Y) / sim.sample_size
            del X
        map(compute_sumstats_for_slice, d.slices())

        alphahat = np.concatenate((alphahat, myalphahat))

    # write output
    def write_output():
        pickle.dump(indices, sim.individuals_file(
                    args.beta_num, args.sample_num, 'wb'), 2)
        pickle.dump(Y, sim.noisy_Y_file(
                    args.beta_num, args.sample_num, 'wb'), 2)
        pickle.dump(alphahat, sim.sumstats_file(
                    args.beta_num, args.sample_num, 'wb'), 2)
    write_output()

def submit(args):
    sim = SumstatSimulation(args.sim_name)
    def submit_beta(beta_num):
        my_args = ['--sim-name', args.sim_name,
                'main',
                '--beta-num', str(beta_num),
                '--sample-num', '$LSB_JOBINDEX']
        outfilepath = \
            sim.path_to_beta(beta_num) + \
            '.sim_sumstats.%I.out'

        bsub.submit(
                ['python', '-u', paths.code + 'sim/sim_sumstats.py'] + my_args,
                outfilepath,
                jobname='simsumstats' + str(beta_num) + '[1-' + str(sim.num_samples_per_beta) + ']',
                # memory_GB=10.5)
                memory_GB=13)
    map(submit_beta, range(1,sim.num_betas + 1))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sim-name', type=str, required=True,
            help='the name of the set of simulation parameters to use from simulations.json')

    main_parser, submit_parser = bsub.add_main_and_submit(parser, main, submit)

    main_parser.add_argument('--beta-num', type=int, required=True,
            help='the 1-based index of the beta to use')
    main_parser.add_argument('--sample-num', type=int, required=True,
            help='the 1-based index specifying which independent sample we\'re on for \
                    a given beta')

    bsub.choose_parser_and_run(parser)
