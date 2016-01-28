from __future__ import print_function, division
import json
import pickle
from pyutils import fs
from dataset import Dataset
import paths


class SumstatSimulation(object):
    def __init__(self, name, path=paths.simulations):
        self.name = name
        self.__dict__.update(
                json.load(open(path + name + '.json')))

    def __str__(self):
        result = ''
        for n in self.__dict__.keys():
            if '__' not in n:
                result += n + '\t' + str(self.__dict__[n])
        return result

    def readable_name(self):
        return '{},{},h2g={},sample_size={}'.format(
                self.dataset,
                self.architecture,
                self.h2g,
                self.sample_size)

    def path_to_refpanel(self):
        return paths.datasets + self.dataset + '/'

    def path(self, create=True):
        path = self.path_to_refpanel() + self.name + '/'
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

    def sumstats_files(self, beta_num):
        for i in range(1, self.num_samples_per_beta+1):
            yield pickle.load(self.sumstats_file(beta_num, i))


if __name__ == '__main__':
    s = SumstatSimulation('test')
    print(s.__dict__)
