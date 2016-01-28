from __future__ import print_function, division
import numpy as np
import pickle, argparse
import matplotlib.pyplot as plt
from experiment import Experiment


def find_truth(estimators):
    return [e for e in estimators if e.method() == 'Truth'][0]

def create_plot(exp, sim, estimators, error_lists):
    mses = ['{:.2e}'.format(np.sum(errors ** 2)) for errors in error_lists]
    plt.figure()
    plt.boxplot(error_lists,
            labels=[
                e.readable_name().replace(',','\n') + '\n' + mse
                for e, mse in zip(estimators,mses)],
            widths=0.75)
    plt.xticks(rotation=35)

    for i, errors in enumerate(error_lists):
        x = np.random.normal(1+i, 0.06, size=len(errors))
        plt.plot(x, errors, 'r.', alpha=0.1)
    plt.axhline(y=0)
    plt.title(sim.readable_name())
    fig = plt.gcf()
    fig.set_size_inches(8, 8)
    fig.subplots_adjust(bottom=0.25)
    fig.savefig(exp.plot_filename(sim), dpi=400)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp_name', type=str, required=True,
            help='the name of the json experiment file whose results we wish to plot')
    args = parser.parse_args()

    exp = Experiment(args.exp_name)

    def plot_results_for(sim):
        truth = find_truth(exp.estimators)
        error_lists = []
        for e in exp.estimators:
            if e == truth: continue
            my_errors = np.empty((0,))
            for beta in range(1, sim.num_betas+1):
                true_results = truth.results(beta, sim)
                my_errors = np.append(my_errors, e.results(beta, sim) - true_results)
            error_lists.append(my_errors)
        create_plot(exp, sim, exp.estimators, error_lists)
    map(plot_results_for, exp.simulations)
