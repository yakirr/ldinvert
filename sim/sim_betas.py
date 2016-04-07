from __future__ import print_function, division
import numpy as np
import pickle, argparse, subprocess
from time import time
from pyutils import bsub
from primitives import SumstatSimulation, Architecture, Dataset
import paths


def main(args):
    np.random.seed(args.beta_num)
    sim = SumstatSimulation(args.sim_name)
    arch = Architecture(sim.architecture)

    beta = np.array([])
    Y = 0
    for chrnum in sim.chromosomes:
        print('chr', chrnum)
        d = Dataset(sim.dataset+'.'+str(chrnum))

        # sample the beta
        mybeta = arch.draw_effect_sizes(sim.dataset, sim.h2g, chrnum, sim.chromosomes)[:,0]
        beta = np.concatenate((beta, mybeta))

        # compute noiseless phenotypes slice by slice
        t0 = time()
        for s in d.slices():
            # X will be N x M
            print(int(time() - t0), ': getting genotypes from file. SNPs', s)
            X = d.get_standardized_genotypes(s)
            print('computing phenotypes. SNPs', s)
            Y += X.dot(mybeta[s[0]:s[1]])
            del X

        # # normalize the Y and the beta to the desired heritability
        # normalization = np.std(Y) / np.sqrt(sim.h2g)
        # if normalization == 0: normalization = 1  # just in case we have some 0s...
        # Y /= normalization
        # mybeta /= normalization

    # write the betas and the noiseless phenotypes
    pickle.dump(beta, sim.beta_file(args.beta_num, 'wb'), 2)
    pickle.dump(Y, sim.noiseless_Y_file(args.beta_num, 'wb'), 2)

def submit(args):
    sim = SumstatSimulation(args.sim_name)
    my_args = ['--sim-name', args.sim_name,
            'main',
            '--beta-num', '$LSB_JOBINDEX']
    outfilepath = sim.path() + \
            '.sim_betas.%I.out'
    bsub.submit(
            ['python', '-u', paths.code + 'sim/sim_betas.py'] + my_args,
            outfilepath,
            jobname='simbetas[1-' + str(sim.num_betas) + ']',
            memory_GB=5)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sim-name', type=str, required=True,
            help='the name of the set of simulation parameters to use from simulations.json')

    main_parser, submit_parser = bsub.add_main_and_submit(parser, main, submit)
    main_parser.add_argument('--beta-num', type=int, required=True,
            help='the 1-based index of the beta to generate')

    bsub.choose_parser_and_run(parser)
