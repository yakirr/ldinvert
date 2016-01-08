from __future__ import print_function, division
import pickle
from meta import Annotation

class LDSC(object):
    def __init__(self, params, dataset_name):
        self.__dict__.update(params)
        self.annotation = Annotation(self.annotation_name, dataset_name)
        self.irs = pickle.load(self.annotation.irsfile())

        self.precompute()

    def precompute(self):
        ldscores_annot = pickle.load(
                self.annotation.ldscores_file_annot(self.Nref,
                    snp_bandwidth=self.bandwidth))
        ldscores_notannot = pickle.load(
                self.annotation.ldscores_file_notannot(self.Nref,
                    snp_bandwidth=self.bandwidth))
        self.ldscores = ldscores_annot + ldscores_notannot
        self.total_ldscore = np.sum(self.ldscores)
        self.avg_ldscore = np.mean(self.ldscores)
        self.X = np.array([
            np.ones(len(ldscores_annot)),
            ldscores_annot,
            ldscores_notannot]).T
        print('total ld score to pathway:', np.sum(self.ldscores_annot))
        print('total ld score to rest of genome:', np.sum(self.ldscores_notannot))

    def run(alphahat, indiv_indices, Y, N):
        chisq = N * alphahat ** 2

        tauhat = (np.mean(chisq) - 1) / (N * self.avg_ldscore)
        w = np.sqrt(np.array([
            1 / (l * (1 + N * tauhat * l)**2) for l in self.ldscores
            ]))
        weighted_X = self.X * w[:,None]
        hat_matrix = np.linalg.solve(weighted_X.T.dot(weighted_X), weighted_X.T)

        weighted_chisq = w * chisq
        solution = hat_matrix.dot(weighted_chisq)
        return solution[1] * len(self.irs) / N




#TODO: merge code below into new LDSC class, then do the same for other methods,
#       then hook up meta/estimator.py to run/run_estimators.py via create_plan.py, which
#           should be renamed to build_analysis.py or plan_experiment.py or something and
#           should be moved to then 'run' folder.
#       consider making 'experiment_name' a global variable on the theory that you should
#           only be doing one experiment at a time and it's a pain in the ass to keep typing
#           its name. What if you want to submit more than one experiment though? maybe not.
