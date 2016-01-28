from __future__ import print_function, division
import json
from primitives import SumstatSimulation
import methods
from pyutils import fs
import paths


class Estimators(object):
    def __init__(self, name):
        self.estimators = []
        for e_def in json.load(open(paths.estimator_sets + name + '.json')):
            method_name = e_def['method']; del e_def['method']
            self.estimators.append(
                    methods.find_method(method_name)(**e_def))

    def __iter__(self):
        return self.estimators.__iter__()


class Simulations(object):
    def __init__(self, name):
        self.simulations = []
        for s_name in json.load(open(paths.simulation_sets + name + '.json')):
            self.simulations.append(SumstatSimulation(s_name))

    def __iter__(self):
        return self.simulations.__iter__()


class Experiment(object):
    def __init__(self, name):
        self.name = name
        self.__dict__.update(json.load(
            open(paths.experiments + name + '.json')))
        self.estimators = Estimators(self.estimator_set)
        self.simulations = Simulations(self.sim_set)

    def results_folder(self, create=True):
        path = paths.results + self.name + '/'
        if create:
            fs.makedir(path)
        return path

    def plot_filename(self, sim):
        return self.results_folder() + sim.readable_name() + '.results.png'

if __name__ == '__main__':
    exp = Experiment('testexp')
    ests = Estimators(exp.estimator_set)
    sims = Simulations(exp.sim_set)

    for s in sims:
        for e in ests:
            print(str(e), str(s.name))
