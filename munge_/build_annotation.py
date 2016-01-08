from __future__ import print_function, division
import numpy as np
import pickle, os, argparse
from pybedtools import BedTool, Interval
from pysnptools.util.intrangeset import IntRangeSet
import pyutils.fs as fs
import genome.utils as gutils
from meta import Annotation, Dataset
import paths

def save_snpbelonging(bedfile, outfile, dataset):
    snps_ucscbed = dataset.snps_bedtool()
    pathway_indicator = snps_ucscbed.intersect(bedfile, c=True) # create indicator vector
    pathway_indices = IntRangeSet(np.flatnonzero(
        np.array(
            [int(snp.name) for snp in pathway_indicator]
        )))
    print('total size of', outfile.name, '=', len(pathway_indices), 'SNPs')
    pickle.dump(pathway_indices, outfile, 2)

def build_annotation(GO_name, window_size_Mb, dataset):
    # Load up pathway
    print('loading pathway')
    pathway = BedTool(paths.GO_dir + GO_name + '.bed').sort()
    annotation = Annotation(GO_name, dataset.name)
    window_size_str = str(window_size_Mb) + 'Mb'

    print(pathway)
    print('total size before window addition:', pathway.total_coverage(), 'bp')

    # Create directory for output
    fs.makedir(annotation.path_to_refdata())

    # Intersect the pathway with the appropriate genome build
    # TODO: this step should be unnecessary if the pathways are correct
    print('removing regions not in reference')
    chromosomes = BedTool(paths.reference + dataset.reference_genome + '.chrlen.bed')
    pathway = chromosomes.intersect(pathway).sort().saveas(
        annotation.ucscbedfilename())

    # compute the flanks
    # TODO: use 1cM instead of 1Mb
    print('computing flanks')
    flanks = pathway.flank(
        genome=dataset.reference_genome,
        b=window_size_Mb*1000000).sort().merge().saveas(
            annotation.flanks_ucscbedfilename(window_size=window_size_str))
    print(flanks)

    # compute the union of the flanks and the pathway
    print('computing union')
    union = pathway.cat(flanks, postmerge=False).sort()
    merged = union.merge().saveas(
        annotation.merged_ucscbedfilename(window_size=window_size_str))
    print(merged)
    print('total size after window addition:', merged.total_coverage(), 'bp')

    # create an indicator vector for whether each SNP in the relevant
    # data set is in the pathway
    print('computing SNP indexset for each region')
    save_snpbelonging(pathway,
            annotation.irsfile(mode='wb'), dataset)
    save_snpbelonging(flanks,
            annotation.flanks_irsfile(mode='wb', window_size=window_size_str), dataset)
    save_snpbelonging(merged,
            annotation.merged_irsfile(mode='wb', window_size=window_size_str), dataset)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--GO_name', type=str, required=True,
            help='The name of the GO file from which to build the annotation')
    parser.add_argument('--dataset_name', type=str, required=True,
            help='The name of the dataset to which to customize the annotation')
    parser.add_argument('--window_size', type=int, required=True,
            help='The amount, in Mb, to add around the pathway to account for LD')
    args = parser.parse_args()
    build_annotation(args.GO_name,
        args.window_size,
        Dataset(args.dataset_name))
