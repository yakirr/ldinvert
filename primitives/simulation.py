from __future__ import print_function, division
import json
from pyutils import fs
from dataset import Dataset
import paths


class SumstatSimulation(object):
    def __init__(self, name, path=paths.simulations):
        self.name = name
        self.__dict__.update(
                json.load(open(path + name + '.json')))

    def long_name(self):
        return self.architecture + ',N=' + str(self.sample_size)

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

    # decided that refpanel size is part of the parameters of each method

if __name__ == '__main__':
    s = SumstatSimulation('test')
    print(s.__dict__)
