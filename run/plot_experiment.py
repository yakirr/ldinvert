from __future__ import print_function, division
import numpy as np
import pickle, argparse
import matplotlib.pyplot as plt
from experiment import Experiment


def create_plot(exp, sim, error_lists, pretty=False):
    print('inside create_plot')
    mses = {e:'{:.2e}'.format(np.mean(error_lists[e]**2)) for e in error_lists.keys()}
    medses = {e:'{:.2e}'.format(np.median(error_lists[e]**2)) for e in error_lists.keys()}

    def label(e):
        if not pretty:
            return e.readable_name().replace(',','\n') + '\n' + mses[e] + '\n' + medses[e]
        else:
            return e.params.pretty_name.replace(',','\n')

    # create box plots
    plt.figure()
    plt.boxplot([error_lists[e].reshape((-1,)) for e in exp.estimators],
        labels=[label(e) for e in exp.estimators],
        sym='',
        widths=0.75)

    # add individual points with jitter
    for i, e in enumerate(exp.estimators):
        # x = np.random.normal(1+i, 0.06, size=len(error_lists[e]))
        x = (1 + i) + np.arange(-0.3, 0.3, 0.6/len(error_lists[e]))
        plt.plot(x, error_lists[e], 'r.', alpha=0.05)
        plt.plot(x, np.mean(error_lists[e], axis=1), 'b.') # to show average for each beta
        plt.plot([1+i], np.mean(error_lists[e]), 'k.')

    # formatting
    plt.xticks(rotation=35)
    plt.axhline(y=0)
    plt.ylabel('error')
    if not pretty:
        plt.title(sim.readable_name())
    fig = plt.gcf()
    fig.set_size_inches(2 + len(error_lists.keys()) * 1.5, 10)
    fig.subplots_adjust(bottom=0.25)

    print('saving figure')
    fig.savefig(exp.plot_filename(sim), dpi=300)
    print('showing figure')
    plt.show()
    fig.gca().set_ylim([-0.1, 0.1])
    print('saving with standardized axes')
    fig.savefig(exp.plot_filename(sim) + '.axes.png', dpi=300)

def write_results(exp, sim, result_lists, truth_lists, error_lists):
    for e in exp.estimators:
        with open(exp.resultstsv_filename(sim, e), 'w') as f:
            print('beta\ttruth\tresult\terror', file=f)
            for beta in range(1, sim.num_betas+1):
                for t,r,err in zip(truth_lists[e][beta-1], result_lists[e][beta-1], error_lists[e][beta-1]):
                    print(beta, t, r, err, sep='\t', file=f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp_name', type=str, required=True,
            help='the name of the json experiment file whose results we wish to plot')
    parser.add_argument('-pretty', action='store_true', default=False,
            help='print short pretty names for methods and leave out mse\'s')
    args = parser.parse_args()

    print('loading experiment', args.exp_name)
    exp = Experiment(args.exp_name)

    if not args.pretty:
        plt.rcParams.update({'axes.titlesize': 'medium'})
        plt.rcParams.update({'font.size': 7})
    else:
        plt.rcParams.update({'axes.titlesize': 'medium'})
        plt.rcParams.update({'font.size': 9})

    print('plotting results')
    def plot_results_for(sim):
        error_lists = {}
        result_lists = {}
        truth_lists = {}
        for e in exp.estimators:
            my_errors = np.empty((sim.num_betas,sim.num_samples_per_beta))
            my_results = np.empty((sim.num_betas,sim.num_samples_per_beta))
            my_truths = np.empty((sim.num_betas,sim.num_samples_per_beta))
            for beta in range(1, sim.num_betas+1):
                my_truths[beta-1] = exp.truth.results(beta, sim)
                my_results[beta-1] = e.results(beta, sim)
                my_errors[beta-1] = my_results[beta-1] - my_truths[beta-1]

            error_lists[e] = my_errors
            result_lists[e] = my_results
            truth_lists[e] = my_truths
        create_plot(exp, sim, error_lists, pretty=args.pretty)
        write_results(exp, sim, result_lists, truth_lists, error_lists)
        with open(exp.purpose_filename(), 'w') as outfile:
            outfile.write(exp.purpose)
    map(plot_results_for, exp.simulations)

