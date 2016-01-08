from __future__ import print_function, division
import numpy as np
import pickle, argparse, subprocess
import pyutils.iter as it
from pyutils import bsub
from meta import Annotation, Dataset
import paths

def main(args):
    dataset = Dataset(args.dataset_name)
    annotation = Annotation(args.set_name, args.dataset_name)
    pathway_indexset = pickle.load(annotation.irsfile())

    # sample random set of 14k individuals (this is done for ldscore and for covariance)
    np.random.seed(0)
    iids = np.random.choice(dataset.N, size=args.Nref, replace=False)

    ldscores_A = np.zeros(dataset.M)
    ldscores_G = np.zeros(dataset.M)
    def compute_ldscores_for_slice(s):
        print(s)
        genotypes = dataset.genotypes_bedfile()[iids, s[0]:s[1]].read()
        genotypes.standardize(); genotypes.standardize()

        # figure out which snps to compute ld scores for in this slice
        indices = (0 if s[0] == 0 else args.snp_bandwidth,
            dataset.M - s[0] if s[1] == dataset.M else
                genotypes.sid_count-args.snp_bandwidth)

        one_over_N = 1 / genotypes.iid_count
        def compute_ldscores_for_snp(m):
            start = max(0, m - args.snp_bandwidth)
            end = min(genotypes.sid_count, m + args.snp_bandwidth)

            # compute the scores
            window = genotypes.val[:, start:end]
            r2_to_snps_in_window = (genotypes.val[:,m].dot(window) / genotypes.iid_count)**2
            r2_to_snps_in_window -= one_over_N

            # add to appropriate ldscore buckets
            pathway_indicator = np.array([s[0] + m_prime in pathway_indexset
                for m_prime in range(start, end)])
            ldscores_A[s[0] + m] += np.sum(r2_to_snps_in_window * pathway_indicator)
            ldscores_G[s[0] + m] += np.sum(
                    r2_to_snps_in_window * np.logical_not(pathway_indicator))

        map(compute_ldscores_for_snp, it.show_progress(range(indices[0], indices[1])))
    map(compute_ldscores_for_slice,
            dataset.slices(buffer_size=args.snp_bandwidth))

    # write output
    with annotation.ldscores_file_annot(args.Nref,
            snp_bandwidth=args.snp_bandwidth,
            mode='wb') as ldscores_A_file:
        pickle.dump(ldscores_A, ldscores_A_file, 2)
    with annotation.ldscores_file_notannot(args.Nref,
            snp_bandwidth=args.snp_bandwidth,
            mode='wb') as ldscores_G_file:
        pickle.dump(ldscores_G, ldscores_G_file, 2)

def submit(args):
    my_args = ['--dataset_name', str(args.dataset_name),
            '--set_name', str(args.set_name),
            '--Nref', str(args.Nref),
            '--snp_bandwidth', str(args.snp_bandwidth),
            'main']
    outfilepath = Annotation(args.set_name, args.dataset_name).path_to_refdata() + \
        '.out.ldscores.' + ','.join(my_args[:-1]).replace('-','')

    cmd = bsub.bsub_command(
            ['python', '-u', paths.code + 'munge/pathway_ldscores.py'] + my_args,
            outfilepath,
            jobname='ldscores,N='+str(args.Nref),
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
    parser.add_argument('--snp_bandwidth', type=int, default=500,
            help='the size of the window to place around the annotation, in SNPs. E.g., 500')

    bsub.add_main_and_submit(parser, main, submit)
    bsub.choose_parser_and_run(parser)
