from __future__ import print_function, division
import json
import pickle
import pandas as pd
from pyutils import fs
from dataset import Dataset
import paths


class SumstatSimulation(object):
    def __init__(self, name, path=paths.simulations):
        self.name = name
        self.__dict__.update(
                json.load(open(path + name + '.json')))
        if self.chromosomes is not None:
            self.__dataset = {
                    c : Dataset(self.dataset+'.'+str(c))
                    for c in self.chromosomes}
            self.by_chrom = True
        else:
            self.__dataset = {
                    0 : Dataset(self.dataset)
                    }
            self.by_chrom = False

    def __str__(self):
        result = ''
        for n in self.__dict__.keys():
            if '__' not in n:
                result += n + '\t' + str(self.__dict__[n])
        return result

    def readable_name(self):
        return '{},{},h2g={},sample_size={},cond={}'.format(
                self.dataset,
                self.architecture,
                self.h2g,
                self.sample_size,
                self.condition_on_covariates)

    def path_to_genotypes(self):
        return self.__dataset.values()[0].path
    def path_to_auxfiles(self):
        if self.by_chrom:
            return self.__dataset.values()[0].auxfiles_path + '../'
        else:
            return self.__dataset.values()[0].auxfiles_path

    def path(self, create=True):
        path = self.path_to_auxfiles() + self.name + '/'
        if create:
            fs.makedir(path)
        return path
    def path_to_beta(self, beta_num, create=True):
        path = self.path(create=create) + 'beta.' + str(beta_num) + '/'
        if create:
            fs.makedir(path)
        return path

    def beta_file(self, beta_num, mode='rb'):
        return open(self.path_to_beta(beta_num) + 'beta', mode)
    def noiseless_Y_file(self, beta_num, mode='rb'):
        return open(self.path_to_beta(beta_num) + 'noiseless_Y.1darray', mode)
    def noisy_Y_file(self, beta_num, index, mode='rb'):
        return open(self.path_to_beta(beta_num) +
                str(index) + '.Y.1darray', mode)
    def individuals_file(self, beta_num, index, mode='rb'):
        return open(self.path_to_beta(beta_num) +
                str(index) + '.indivs_indices.1darray', mode)
    def sumstats_file(self, beta_num, index, mode='rb'):
        return open(self.path_to_beta(beta_num) +
                str(index) + '.alphahat', mode)

    def sumstats_aligned_to_refpanel(self, beta_num, refpanel, chrnum):
        to_flip = self.__dataset[chrnum].snp_consistency_vector(refpanel)

        for alphahat in self.sumstats_files(beta_num):
            alphahat[to_flip] *= -1
            yield alphahat

    def sumstats_files(self, beta_num):
        for i in range(1, self.num_samples_per_beta+1):
            yield pickle.load(self.sumstats_file(beta_num, i))


if __name__ == '__main__':
    s = SumstatSimulation('test')
    print(s.__dict__)
