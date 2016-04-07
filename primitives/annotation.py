from __future__ import print_function, division
import numpy as np
import pandas as pd
import paths

class Annotation(object):
    def __init__(self, stem, chrnum=None):
        self.stem = stem
        if chrnum is None:
            self.split_by_chr = False
        else:
            self.split_by_chr = True
            self.chrnum = chrnum

    @property
    def filestem(self):
        if self.split_by_chr:
            return '{}{}.{}'.format(paths.annotations, self.stem, self.chrnum)
        else:
            return '{}{}'.format(paths.annotations, self.stem)

    @property
    def filename(self):
        return self.filestem + '.annot.gz'

    @property
    def norm_filename(self):
        return self.filestem + '.norm'

    def get_vectors(self, refpanel):
        df = pd.read_csv(self.filename, compression='gzip', header=0,
                sep='\t')
        names = df.columns.values[4:]
        refpanel_bim = pd.read_csv(refpanel.genotypes_bedfile.filename + '.bim',
                sep='\t',
                names=['CHR', 'SNP', 'cM', 'BP', 'A1', 'A2'])
        result = refpanel_bim.merge(df, how='left', on=['SNP'])
        return {
                n: result[n].values
                for n in names}
    def get_vector(self, name):
        return self.get_vectors()[name]

    def get_norms(self):
        df = pd.read_csv(self.norm_filename, header=0, sep='\t')
        return {
                n : df[n][0]
                for n in df.columns.values}
    def get_norm(self, name):
        return self.get_norms()[name]


if __name__ == '__main__':
    from dataset import Dataset
    a = Annotation('/groups/price/yakir/data/signed_annot/maf5p')
    refpanel = Dataset('UK10Khg19.22')
    print(a.get_vectors(refpanel, 22))
