from __future__ import print_function, division
import numpy as np
import json
from genome import GenomicSubset, SnpSubset
from dataset import Dataset
from annotation import Annotation
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
    def __init__(self, name):
        self.name = name
        self.__dict__.update(
                json.load(open(paths.architectures + name + '.json')))

    # returns an array that is genome size x num_samples
    def draw_effect_sizes(self, dataset_name, Eh2g, chrnum, allchroms, num_samples=1):
        d = Dataset(dataset_name + '.' + str(chrnum))

        # load annotations
        annotations = {}
        for annot_file in self.annot_files:
            a = Annotation(annot_file, chrnum)
            annotations.update(a.get_vectors(d))

        # figure out total norms of the annotations
        norms = {}
        for annot_file in self.annot_files:
            for c in allchroms[:1]:
                a = Annotation(annot_file, c)
                norms.update(a.get_norms())
            for c in allchroms[1:]:
                a = Annotation(annot_file, c)
                for n, norm in a.get_norms().items():
                    norms[n] += norm
        print(annotations)
        print(norms)

        total_variance = np.sum([norms[n] * mu**2 for n, mu in self.mean_effects.items()])
        total_variance += np.sum([
            norms[n] * e['sigma2'] for n, e in self.variance_effects.items()])
        normalization = Eh2g / total_variance

        class local:
            result = np.zeros((d.M, num_samples))
        def add_variance_effect((n, e)):
            per_snp_variance = e['sigma2'] * normalization
            print('norm of', n, 'is', norms[n])
            print('per-SNP variance is', per_snp_variance)
            print('contributing heritability of', per_snp_variance * norms[n], 'over all chr')
            esd = EffectSizeDist(e['effectdist'])
            local.result[np.flatnonzero(annotations[n])] += esd.draw_effect_sizes(
                    np.count_nonzero(annotations[n]), num_samples, per_snp_variance)
        map(add_variance_effect, self.variance_effects.items())

        def add_mean_effect((n, mu)):
            real_mu = mu * np.sqrt(normalization)
            print('norm of', n, 'is', norms[n])
            print('per-SNP mu is', real_mu)
            print('contributing heritability of', real_mu**2 * norms[n], 'over all chr')
            local.result += (annotations[n] * real_mu)[:,None]
        map(add_mean_effect, self.mean_effects.items())

        return local.result

if __name__ == '__main__':
    np.random.seed(0)
    esd = EffectSizeDist('sparse')
    print('below should be a sparse 5x10 matrix with large non-zero values')
    print(esd.draw_effect_sizes(5, 10, 100))
    print()

    a = Architecture('test')
    betas = a.draw_effect_sizes('UK10Khg19', 0.5, 22, [22], num_samples=2)
    print(betas[290:310])
