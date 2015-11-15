from __future__ import print_function, division
import hyperparams as hp
import pyutils
import argparse
import os, sys, subprocess
import estimators
import pickle

parser = argparse.ArgumentParser()

def output_filename(annotation, method):
    return method + ':' + annotation

def analyze(method, first_beta=1):
    def get_truth(beta_num):
        return estimators.truth(
            pickle.load(hp.beta_file(beta_num)))

    def analyze_beta(beta_num, method):
        def analyze_sample(beta_num, index, method):
            alphahat = pickle.load(hp.sumstats_file(beta_num, index))
            indiv_indices = pickle.load(hp.individuals_file(beta_num, index))
            Y = pickle.load(hp.noisy_Y_file(beta_num, index))
            N = len(Y)
            return estimators.methods[method](alphahat, indiv_indices, Y, N)

        index = 1
        results = []
        while True:
            try:
                test = hp.noisy_Y_file(beta_num, index); test.close()
            except IOError:
                break
            results.append(analyze_sample(beta_num, index, method))
            index += 1
            print('.', end='.')
            sys.stdout.flush()
        print()
        return results

    beta_num = first_beta
    beta_nums, truths, results = [], [], []
    while os.path.exists(hp.path_to_samplesize_dir(beta_num, create=False)):
        print(beta_num)
        truth = get_truth(beta_num)
        print(truth)
        these_results = analyze_beta(beta_num, method)
        results.extend(these_results)
        truths.extend([truth] * len(these_results))
        beta_nums.extend([beta_num] * len(these_results))
        beta_num += 1
    return beta_nums, truths, results

def write_results(beta_nums, truths, results):
    outfile = open(hp.path_to_results_dir() + 'results:' + \
            output_filename(args.annotation_name, args.method_name), 'w')
    outfile.write('beta_num\ttruth\tresult\n')
    for beta_num, truth, result in zip(beta_nums, truths, results):
        outfile.write('{}\t{}\t{}\n'.format(beta_num, truth, result))
    outfile.close()

def main(args):
    hp.load()
    estimators.initializers[args.method_name](args.annotation_name)
    write_results(*analyze(args.method_name))

def submit(args):
    hp.load(printall=False)
    def submit_method(method_name):
        my_args = ['--annotation_name', args.annotation_name,
                'main',
                '--method_name', method_name]
        outfilepath = hp.path_to_results_dir() + '.out.assess:' + output_filename(
                args.annotation_name, method_name)

        cmd = pyutils.configs.bsub_command(
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
    parser.add_argument('--annotation_name', type=str, required=True,
            help='the name of the annotation in which to assess enrichment')

    subparsers = parser.add_subparsers()

    subparser_main = subparsers.add_parser('main')
    subparser_main.add_argument('--method_name', type=str, required=True,
            help='the name of the estimator to use. names specified in est/estimators.py')
    subparser_main.set_defaults(_func=main)

    subparser_submit = subparsers.add_parser('submit')
    subparser_submit.set_defaults(_func=submit)

    args, _ = parser.parse_known_args()
    pyutils.configs.print_vars(args)
    print()

    args._func(args)
