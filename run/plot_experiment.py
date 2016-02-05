from __future__ import print_function, division
import numpy as np
import pickle, argparse
import matplotlib.pyplot as plt
from experiment import Experiment


def create_plot(exp, sim, error_lists):
    mses = {e:'{:.2e}'.format(np.mean(error_lists[e]**2)) for e in error_lists.keys()}
    def label(e):
        return e.readable_name().replace(',','\n') + '\n' + mses[e]

    sorted_estimators = sorted(error_lists.keys(), key=lambda e: e.readable_name())

    # create box plots
    plt.figure()
    plt.boxplot([error_lists[e] for e in sorted_estimators],
        labels=[label(e) for e in sorted_estimators],
        widths=0.75)

    # add individual points with jitter
    for i, e in enumerate(sorted_estimators):
        x = np.random.normal(1+i, 0.06, size=len(error_lists[e]))
        plt.plot(x, error_lists[e], 'r.', alpha=0.1)

    # formatting
    plt.xticks(rotation=35)
    plt.axhline(y=0)
    plt.title(sim.readable_name())
    fig = plt.gcf()
    fig.set_size_inches(len(error_lists.keys()) * 1.5, 8)
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
        error_lists = {}
        for e in exp.estimators:
            my_errors = np.empty((0,))
            for beta in range(1, sim.num_betas+1):
                true_results = exp.truth.results(beta, sim)
                my_errors = np.append(my_errors, e.results(beta, sim) - true_results)
            error_lists[e] = my_errors
        create_plot(exp, sim, error_lists)
        with open(exp.purpose_filename(), 'w') as outfile:
            outfile.write(exp.purpose)
    map(plot_results_for, exp.simulations)
