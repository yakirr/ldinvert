from __future__ import print_function, division
import argparse, pickle, subprocess
import os, sys
from pyutils import bsub
from meta import Annotation, Dataset, Estimator, Simulation

def analyze(estimator, dataset, simulation, first_beta=1):
    def analyze_beta(beta_num):
        def analyze_sample(beta_num, index):
            alphahat = pickle.load(simulation.sumstats_file(dataset.name, beta_num, index))
            indiv_indices = pickle.load(
                    simulation.individuals_file(dataset.name, beta_num, index))
            Y = pickle.load(simulation.noisy_Y_file(dataset.name, beta_num, index))
            return estimator.run(alphahat, indiv_indices, Y, simulation.N)

        index = 1
        results = []
        while True:
            try:
                test = simulation.noisy_Y_file(dataset.name, beta_num, index); test.close()
            except IOError:
                break
            results.append(analyze_sample(beta_num, index))
            index += 1
            print('.', end='.')
            sys.stdout.flush()
        print()
        return results

    beta_num = first_beta
    beta_nums, truths, results = [], [], []
    while os.path.exists(
            simulation.path_to_samplesize_dir(dataset.name, beta_num, create=False)):
        print(beta_num)
        beta = pickle.load(simulation.beta_file(dataset.name, beta_num))
        truth = estimator.truth(beta)
        print(truth)
        these_results = analyze_beta(beta_num)
        results.extend(these_results)
        truths.extend([truth] * len(these_results))
        beta_nums.extend([beta_num] * len(these_results))
        beta_num += 1
    return beta_nums, truths, results

def write_results(estimator, beta_nums, truths, results):
    with estimator.results_file(mode='w') as outfile:
        outfile.write('beta_num\ttruth\tresult\n')
        for beta_num, truth, result in zip(beta_nums, truths, results):
            outfile.write('{}\t{}\t{}\n'.format(beta_num, truth, result))

def main(args):
    print('main')
    d = Dataset(args.dataset_name)
    s = Simulation(args.simulation_name)
    e = Estimator(args.estimator_name, experiment_name, dataset.name)
    # write_results(*analyze(e, d, s))

def submit(args):
    hp.load(printall=False)
    def submit_method(method_name):
        my_args = ['--annotation_name', args.annotation_name,
                'main',
                '--method_name', method_name]
        outfilepath = hp.path_to_results_dir() + '.out.assess:' + output_filename(
                args.annotation_name, method_name)

        cmd = bsub.bsub_command(
                ['python', '-u', hp.paths.code + 'est/parallel_assess.py'] + \
                        my_args + hp.to_command_line(),
                outfilepath,
                jobname='assess_' + method_name,
                memory_GB=8)
        print(' '.join(cmd))
        print(outfilepath)
        subprocess.call(cmd)
    map(submit_method, estimators.methods.keys())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset_name', type=str, required=True,
            help='the name of the dataset whose sumstats should be analyzed')
    parser.add_argument('--simulation_name', type=str, required=True,
            help='the name of the simulation parameter set to be analyzed')
    parser.add_argument('--estimator_name', type=str, required=True,
            help='the full name of the estimator to be used')

    args, _ = parser.parse_known_args()
    bsub.print_namespace(args)
    print()

    main(args)
