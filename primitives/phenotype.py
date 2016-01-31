from __future__ import print_function, division
import numpy as np
import json
from genome import GenomicSubset, SnpSubset
from dataset import Dataset
import paths


class EffectSizeDist(object):
    def __init__(self, name):
        self.name = name
        self.__dict__.update(json.load(
            open(paths.architectures + 'effect_distributions/' + name + '.json')))

    # returns an array that is num_snps x num_samples
    def draw_effect_sizes(self, num_snps, num_samples, per_snp_variance):
        if self.causal_dist == 'normal':
            result = np.zeros((num_snps, num_samples))
            causal_mask = np.random.multinomial(
                    1,
                    [self.p_causal, 1-self.p_causal],
                    size=(num_samples,num_snps))[:,:,0].T.astype(bool)
            num_causal = np.sum(causal_mask)
            result[causal_mask] = \
                    np.sqrt(per_snp_variance / self.p_causal) * np.random.randn(num_causal)
            return result
        else:
            print('error, unsupported causal effect-size distribution')
            return None

class Architecture(object):
    def __init__(self, name, categories=None):
        self.name = name
        if categories:
            self.categories = categories
        else:
            self.categories = json.load(
                    open(paths.architectures + name + '.architecture.json'))

    # returns an array that is genome size x num_samples
    def draw_effect_sizes(self, dataset_name, Eh2g, num_samples=1):
        d = Dataset(dataset_name)
        total_variance = float(sum([i['total_variance'] for i in self.categories.values()]))
        result = np.zeros((d.M, num_samples))

        def add_category_effects((gs, i)):
            ss = SnpSubset(gs, d)
            per_snp_variance = \
                    float(Eh2g * i['total_variance']) / (total_variance * ss.num_snps())
            print('size of', gs.name, 'is', ss.num_snps())
            print('per-SNP variance is', per_snp_variance)
            esd = EffectSizeDist(i['effectdist'])
            result[ss.irs] += esd.draw_effect_sizes(
                    ss.num_snps(), num_samples, per_snp_variance)
        map(add_category_effects, [(GenomicSubset(c), i) for c, i in self.categories.items()])

        return result

if __name__ == '__main__':
    np.random.seed(0)
    esd = EffectSizeDist('sparse')
    print('below should be a sparse 5x10 matrix with large non-zero values')
    print(esd.draw_effect_sizes(5, 10, 100))
    print()

    a = Architecture('tiny_test')
    betas = a.draw_effect_sizes('tinyGERA', 0.5, num_samples=2)
    print()
    print('below should be a 20x2 matrix whose first 10 rows are mostly 0 or small numbers \
        and whose second 10 rows are all non-zero numbers that are larger')
    print(betas[290:310])
