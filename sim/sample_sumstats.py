from __future__ import print_function, division
import numpy as np
import pickle
import argparse
import subprocess
import pyutils.configs
import genome.utils as gutils
import pyutils.iter as it
import sim.utils as sutils
import hyperparams as hp

parser = argparse.ArgumentParser()

def main(args):
    np.random.seed(args.beta_num)
    hp.load()

    # read in noiseless phenotypes
    Y = pickle.load(hp.noiseless_Y_file(args.beta_num))

    # choose individuals and create ensemble of Ys
    indices = np.random.choice(Y.shape[0], size=(args.num_samples, hp.sumstats_params.N))
    Ys = Y[indices]

    # compute how much noise to add
    sigma2e = 1 - (hp.beta_params.h2gA + hp.beta_params.h2gG)
    print('adding noise. sigma2e =', sigma2e)
    Ys += np.sqrt(sigma2e) * np.random.randn(*Ys.shape)

    chrnum_to_alphahats = {}
    def compute_sumstats_for_chr(chrnum):
        print('computing sumstats. chrom', chrnum)
        # Xs will be N x M
        print('\tgetting genotypes from file')
        Xs = sutils.get_standardized_genotypes(chrnum)

        print('\tcomputing sumstats for each simulation')
        chrnum_to_alphahats[chrnum] = np.array([
                Xs[indices[sample_num]].T.dot(Ys[sample_num])
                for sample_num in range(args.num_samples)
                ]) / hp.sumstats_params.N
    map(compute_sumstats_for_chr, hp.chromosomes())

    # write output
    def write_output_for_samplenum(sample_num):
        pickle.dump(indices[sample_num],
                hp.individuals_file(args.beta_num, sample_num + 1, 'wb'), 2)
        pickle.dump(Ys[sample_num],
                hp.noisy_Y_file(args.beta_num, sample_num + 1, 'wb'), 2)
        alphahat_i = {
                chrnum:alphahats[sample_num]
                for chrnum, alphahats in chrnum_to_alphahats.items()}
        pickle.dump(alphahat_i,
                hp.sumstats_file(args.beta_num, sample_num + 1, 'wb'), 2)
    map(write_output_for_samplenum, range(args.num_samples))

def submit(args):
    hp.load(printall=False)
    my_args = ['main',
            '--beta_num', '$LSB_JOBINDEX',
            '--num_samples', str(100)]
    outfilepath = hp.path_to_results_dir() + 'beta.%I/out'

    cmd = pyutils.configs.bsub_command(
            ['python', '-u', hp.paths.code + 'sim/sample_sumstats.py'] + \
                    my_args + hp.to_command_line(),
            outfilepath,
            jobname='samplesumstats[1-20]')
    print(' '.join(cmd))
    print(outfilepath)
    subprocess.call(cmd)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    subparser_main = subparsers.add_parser('main')
    subparser_main.add_argument('--beta_num', type=int, required=True,
            help='the index specifying which beta to use')
    subparser_main.add_argument('--num_samples', type=int, required=True,
            help='the number of random paris of (individuals, noise) to generate')
    subparser_main.set_defaults(_func=main)

    subparser_submit = subparsers.add_parser('submit')
    subparser_submit.set_defaults(_func=submit)

    args, _ = parser.parse_known_args()
    pyutils.configs.print_vars(args)
    print()

    args._func(args)
