from __future__ import print_function, division
import numpy as np
import glob
import hyperparams as hp
import matplotlib.pyplot as plt


def results_text_filename():
    return 'biases:' + hp.pathway.name + '.' + hp.dataset.name + ':' +\
        hp.results_dirname() + '.tsv'

def results_plot_filename():
    return 'biases:' + hp.pathway.name + '.' + hp.dataset.name + ':' +\
        hp.results_dirname() + '.png'

def parse_file(path_to_file):
    def parse(line):
        fields = line.split('\t')
        return np.array([int(fields[0]), float(fields[1]), float(fields[2])])

    with open(path_to_file) as f:
        f.readline()
        [betas, truths, guesses] = np.array([parse(line) for line in f]).T
    return betas, truths, guesses

def process_results_file(path_to_file):
    directory = '/'.join(path_to_file.split('/')[:-1])
    method_name = path_to_file.split('/')[-1].split(':')[1]
    annotation = path_to_file.split('/')[-1].split(':')[2]
    betas, truths, guesses = parse_file(path_to_file)

    uncond_mse = np.mean((guesses - truths)**2) # = expectation over beta of MSE_beta
    uncond_bias = np.mean(guesses - truths) # = expectation over beta of bias_beta
    uncond_var = np.var(guesses) # != expectation over beta of var_beta
    Evar = np.mean([np.var(guesses[betas == beta])
        for beta in betas])
    truthvar = np.var(truths)

    return (method_name, annotation, uncond_mse, uncond_bias,
            uncond_var, Evar, truthvar, guesses - truths)

def create_plot(results):
    plt.figure()
    plt.boxplot([result[7] for result in results],
            labels=[result[0] + '\n' +
                '{:.3e}'.format(result[2]) + '\n' +
                result[1] for result in results],
            widths=0.75)
    plt.xticks(rotation=35)

    for i, result in enumerate(results):
        y = result[7]
        x = np.random.normal(1+i, 0.06, size=len(y))
        plt.plot(x, y, 'r.', alpha=0.1)
    plt.axhline(y=0)
    plt.title('enriched: ' + hp.pathway.name + '.' + hp.dataset.name + '\n' +
            'h2g=' + str(hp.beta_params.h2gG + hp.beta_params.h2gA) + ', '
            'h2gA=' + str(hp.beta_params.h2gA) + ', ' +
            'p_causal(G)=' + str(hp.beta_params.pG) + ', ' +
            'p_causal(A)=' + str(hp.beta_params.pA)
            )
    plt.gcf().set_size_inches(8, 8)
    plt.gcf().subplots_adjust(bottom=0.25)
    plt.savefig(hp.paths.aggregate_results + results_plot_filename(), dpi=400)

if __name__ == '__main__':
    hp.load()
    results = map(process_results_file,
            sorted(glob.glob(hp.path_to_results_dir() + 'results:*')))

    with open(hp.paths.aggregate_results + results_text_filename(), 'w') as f:
        f.write('method\tannotation assessed\tunconditional MSE\n')
        def write_result():
            (method_name, annotation, uncond_mse,
                    uncond_bias, uncond_var, Evar, truthvar, biases) = result
            f.write(method_name + '\t' + str(annotation) + '\t'+ str(uncond_mse))
            f.write('\n')
        map(write_result, results)

    create_plot(results)
