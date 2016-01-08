from __future__ import print_function, division
import numpy as np
from pysnptools.snpreader import Bed

genotypes_path = '/groups/price/hilary/data/'
input_dataset = 'WT1QC'
output_dataset = 'WT1QC_fewnans'

for chrom in range(1, 23):
    print(chrom)
    genotypes_on_disk = Bed(genotypes_path + input_dataset + '/all.' + str(chrom))
    print(genotypes_on_disk)
    genotypes = genotypes_on_disk.read()

    print('removing SNPs with more than 20% missingness...')
    mask = ~(np.sum(np.isnan(genotypes.val), axis=0) / genotypes.val.shape[0] > 0.2)
    print(np.sum(mask), 'out of', genotypes.sid_count, 'remaining')
    del genotypes

    new_genotypes_on_disk = genotypes_on_disk[:, mask]
    new_genotypes = new_genotypes_on_disk.read()
    Bed.write(genotypes_path + output_dataset + '/all.' + str(chrom), new_genotypes)
