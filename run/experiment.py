from __future__ import print_function, division
import itertools
import json
from primitives import SumstatSimulation
import methods
from pyutils import fs
import paths


class Estimators(object):
    def __init__(self, estimators_json):
        self.estimators = []
        for e_def in estimators_json:
            self.estimators.append(Estimators.create_estimator_from_json(e_def))

    def __iter__(self):
        return self.estimators.__iter__()

    @classmethod
    def create_estimator_from_json(cls, json_entry):
        method_name = json_entry['method']; del json_entry['method']
        print(method_name)
        return methods.find_method(method_name)(**json_entry)


class Simulations(object):
    def __init__(self, simulation_names):
        self.simulations = []
        for sim_name in simulation_names:
            self.simulations.append(SumstatSimulation(sim_name))

    def __iter__(self):
        return self.simulations.__iter__()


class Experiment(object):
    def __init__(self, name):
        self.name = name
        self.__dict__.update(json.load(
            open(paths.experiments + name + '.json')))
        self.estimators = Estimators(self.estimators)
        self.truth = Estimators.create_estimator_from_json(self.truth)
        self.simulations = Simulations(self.simulations)

    def results_folder(self, create=True):
        path = paths.results + self.name + '/'
        if create:
            fs.makedir(path)
        return path

    def plot_filename(self, sim):
        return self.results_folder() + sim.readable_name() + '.results.png'

    def resultstsv_filename(self, sim, est):
        return self.results_folder() + sim.readable_name() + '.' + \
                est.readable_name() + '.results.tsv'

    def purpose_filename(self):
        return self.results_folder() + 'purpose.txt'

    def estimators_and_truth(self):
        return itertools.chain(self.estimators, [self.truth])

if __name__ == '__main__':
    exp = Experiment('2016.02.04_ldinvert_Rbiasadjustment')

    for s in exp.simulations:
        for e in exp.estimators_and_truth():
            print(str(e), str(s.name))
