from __future__ import print_function, division
import numpy as np
import pandas as pd
from scipy import stats

print('reading')
path = '/groups/price/yakir/data/sumstats/'
refpanel_bim = '/groups/price/ldsc/reference_files/UK10K/plink_files/ALSPAC.TWINSUK.QC.1.bim'
ss_pos = pd.read_csv(path + 'CD.sumstats.bim', sep='\t',
        names=['chr', 'rsid', 'Morgans', 'bp', 'A1', 'A2'],
        header=False)

refpos = pd.read_csv(refpanel_bim, sep='\t',
        names=['chr', 'rsid', 'Morgans', 'bp', 'A1', 'A2'],
        header=False)

merged = ss_pos.merge(refpos, how='left', on=['rsid'])
ss_pos_ind = np.logical_not(np.isnan(merged['chr_y']))
ss_pos_flip = (merged['A1_x'] == merged['A1_y']) & \
        (merged['A2_x'] == merged['A2_y'])

# 1. load ldblocks bedtool
# 2. load A
# 3a. load sumstats bim
# 3b. load sumstats file
# 4. load refpanel bim
# 5. compute refpanel mask
# 6. compute sumstats mask
# 7. for each ldblock overlapping A:
#   compute irss for ldblock in sumstats and in refpanel
#   let X = refpanel[refpanel_mask & refpanel_ldblock_irs]
#   let alphahat = sumstats[sumstats_mask & sumstats_ldblock_irs]
#   compute R
#   compute Rri
#   etc.


