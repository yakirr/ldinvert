from __future__ import print_function, division
import numpy as np
import pickle, argparse, subprocess
import sparse.blockdiag as bd
from pyutils import bsub
from meta import Annotation, Dataset
import paths


def main(args):
    annotation = Annotation(args.set_name, args.dataset_name)
    dataset = Dataset(args.dataset_name)
    merged_intrangeset = pickle.load(annotation.merged_irsfile(args.window_size))

    # sample random set of 14k individuals (this is done for ldscore and for covariance)
    np.random.seed(0)
    iids = np.random.choice(dataset.N, size=args.Nref, replace=False)

    covariance_matrices = {}
    for r in merged_intrangeset.ranges():
        print('range', r)

        range_data_on_disk = dataset.genotypes_bedfile()
        range_data = range_data_on_disk[iids, r[0]:r[1]].read()
        print('read in all snps in region. shape =', range_data.val.shape)

        print('standardizing...')
        # we need the second call because the first one standardized the variance
        # before substituting 0s for the nans.
        range_data.standardize(); range_data.standardize()

        print('computing covariance...')
        covariance_matrices[r] = \
                range_data.val.T.dot(range_data.val) / range_data.iid_count

    print('assembling into BlockDiagArray...')
    bda = bd.BlockDiagArray(covariance_matrices)

    print('saving...')
    pickle.dump(bda, annotation.covariance_around_annot_file(
        args.Nref,
        banded=False,
        window_size=args.window_size,
        mode='wb'), 2)

def submit(args):
    my_args = ['--dataset_name', str(args.dataset_name),
            '--set_name', str(args.set_name),
            '--Nref', str(args.Nref),
            '--window_size', args.window_size] + \
            (['--banded'] if args.banded else []) + \
            ['main']
    outfilepath = Annotation(args.set_name, args.dataset_name).path_to_refdata() + \
        '.out.covariance.' + ','.join(my_args[:-1]).replace('-','')

    cmd = bsub.bsub_command(
        ['python', '-u', paths.code + 'munge/pathway_covariance.py'] + my_args,
        outfilepath,
        memory_GB=8)
    print(' '.join(cmd))
    print(outfilepath)
    subprocess.call(cmd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--set_name', type=str, required=True,
            help='the name of the subset of the genome to analyze. e.g., "50" or "99".')
    parser.add_argument('--dataset_name', type=str, required=True,
            help='the name of the dataset of genotypes to use')
    parser.add_argument('--Nref', type=int, required=True,
            help='the sample size to use for the estimation')
    parser.add_argument('--window_size', type=str, default='1Mb',
            help='the size of the window to place around the annotation. E.g., "1Mb"')
    parser.add_argument('--banded', action='store_true',
            help='whether to band the estimated covariance matrix. ' +
            '0 means no banding, and anything else means yes banding.')

    bsub.add_main_and_submit(parser, main, submit)
    bsub.choose_parser_and_run(parser)
