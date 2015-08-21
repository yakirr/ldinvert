from __future__ import print_function, division
import os
import numpy as np
import pickle
from pybedtools import BedTool, Interval
from pysnptools.util.intrangeset import IntRangeSet
import pyutils.fs as fs
import partition.paths as pathto
import genome.utils as gutils

dataset_name, genome_build = 'WT1QC_fewnans', 'hg18'
pathway_name = '99'
output_path = pathto.pathways_with_flanks + pathway_name + '.' + dataset_name
window_size_Mb = 500000
chr22only = False

def snp_belonging(region):
    chrom_num = region.chrom[3:]
    snps_in_chrom = BedTool(pathto.genotypes + dataset_name + '/all.' + chrom_num + '.ucscbed')
    indicator = snps_in_chrom.intersect(BedTool([region]), c=True)
    return IntRangeSet(np.flatnonzero(
            np.array(
                [int(snp.name) for snp in indicator]
            )))

def save_snpbelonging(bedfile, name):
    result = {}
    for region in bedfile:
        new_region = Interval(region.chrom, region.start, region.end)
        result[new_region] = snp_belonging(new_region)
    print('total size of', name, '=', gutils.total_size_snps(result.values()), 'SNPs')
    pickle.dump(result, open(output_path + '/' + name + '.regions_to_indexsets', 'wb'), 2)

# Load up pathway
print('loading pathway')
pathway = BedTool('/groups/price/yakir/data/GO/' + pathway_name + '.bed').sort()

# These two lines causes the script to restrict the pathway to chromosome 22
if chr22only:
    pathway = pathway.filter(lambda x: x.chrom == 'chr22').saveas()

print(pathway)
print('total size before window addition:', pathway.total_coverage(), 'bp')

# Create directory for output
fs.makedir(output_path)

# Intersect the pathway with the appropriate genome build
# TODO: this step should be unnecessary if the pathways are correct
print('removing regions not in reference')
chromosomes = BedTool(pathto.reference + genome_build + '.chrlen.bed')
pathway = chromosomes.intersect(pathway).sort().saveas(output_path + '/pathway.bed')

# compute the flanks
# TODO: use 1cM instead of 1Mb
print('computing flanks')
flanks = pathway.flank(genome=genome_build, b=window_size_Mb).sort().merge().saveas(output_path + '/flanks.bed')
print(flanks)

# compute the union of the flanks and the pathway
print('computing union')
union = pathway.cat(flanks, postmerge=False).sort()
merged = union.merge().saveas(output_path + '/merged.bed')
print(merged)
print('total size after window addition:', merged.total_coverage(), 'bp')

# create an indicator vector for whether each SNP in the relevant
# data set is in the pathway
print('computing SNP indexset for each region')
save_snpbelonging(pathway, 'pathway')
save_snpbelonging(flanks, 'flanks')
save_snpbelonging(merged, 'merged')
