
import config

import operator
import os
import itertools

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from faststructure_parse_results import parse_faststructure_results
from passport import get_sample_passports
from pop_building import get_classifications_for_classification_key_path


def plot_marginal_likelihoods(results, plot_path):
    marginal_likelihoods = [(k, res['marginal_likelihood']) for k, res in results.items()]
    marginal_likelihoods = sorted(marginal_likelihoods, key=operator.itemgetter(1))
    ks, marginal_likelihoods = zip(*marginal_likelihoods)

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)
    axes.scatter(ks, marginal_likelihoods)
    axes.set_xlabel('K')
    axes.set_ylabel('Marginal likelihood')

    fig.tight_layout()
    fig.savefig(str(plot_path))


def _plot_admixture_compositions_per_sample(admixtures, plot_path,
                                            pops_for_samples=None,
                                            pop_order=None):
    _plot_admixtures(admixtures, plot_path,
                     pops_for_samples=pops_for_samples, pop_order=pop_order)


def _calc_admixtures_per_pop(admixtures, pops_for_samples):
    pop_col = [pops_for_samples[sample] for sample in admixtures.index]
    admixtures_per_sample = admixtures.copy()
    admixtures_per_sample['pop'] = pop_col

    admixtures_per_pop = admixtures_per_sample.groupby('pop').mean()

    assert numpy.allclose(admixtures_per_pop.sum(axis=1), 1)
    return admixtures_per_pop


def plot_admixtures(admixtures, plot_path=None, axes=None, sample_order=None,
                    pops_for_samples=None, pop_order=None,
                    composition_classes_colors=None):

    if composition_classes_colors is None:
        composition_classes_colors = {}

    if sample_order is None:
        for ancestral_pop_idx in range(admixtures.shape[1]):
            values = admixtures.iloc[:, ancestral_pop_idx]
            sorted_samples = sorted(admixtures.index, key=lambda sample: values.loc[sample])
            admixtures = admixtures.loc[sorted_samples, :]
    else:
        sorted_samples = sorted(admixtures.index, key=lambda sample: sample_order.index(sample))
        admixtures = admixtures.loc[sorted_samples, :]

    if pops_for_samples:
        if pop_order is None:
            pop_order = sorted(set(pops_for_samples.values()))

        pop_order = {pop: idx for idx, pop in enumerate(pop_order)}

        pop_sort_key=lambda x: pop_order[pops_for_samples.get(x)]
        sorted_samples = sorted(admixtures.index, key=pop_sort_key)
        admixtures = admixtures.loc[sorted_samples, :]

    samples = list(admixtures.index)
    bar_edges = numpy.arange(len(samples) + 1)
    bar_centers = (bar_edges[:-1] + bar_edges[1:]) / 2
    width = (bar_edges[1] - bar_edges[0]) * 0.9

    axes.set_facecolor('white')

    if axes is None:
        fig_width = 0.2 * admixtures.shape[0]
        fig = Figure((fig_width, 4))
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)
        savefig = True
    else:
        savefig = False

    ancestral_pops = list(admixtures.columns)
    bottom = None
    for ancestral_pop_idx in range(admixtures.shape[1]):
        ancestral_pop = ancestral_pops[ancestral_pop_idx]
        height = admixtures.iloc[:, ancestral_pop_idx]

        color = composition_classes_colors.get(ancestral_pop, None)

        axes.bar(bar_centers, height, width=width, bottom=bottom, color=color, label=ancestral_pop)

        if bottom is None:
            bottom = height
        else:
            bottom += height

    if pops_for_samples:
        y_min, y_max = axes.get_ylim()
        pops = [pops_for_samples[sample] for sample in admixtures.index]
        start_idx = 0
        for pop, grouped_pops in itertools.groupby(pops):
            n_samples = len(list(grouped_pops))
            end_idx = start_idx + n_samples

            start_pos = bar_edges[start_idx]
            end_pos = bar_edges[end_idx]
            middle_pos = (start_pos + end_pos) / 2
            axes.vlines([start_pos, end_pos], ymin=y_min, ymax=y_max)
            axes.text(middle_pos, y_max, str(pop), rotation=45)

            start_idx = end_idx

    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)

    xtick_labels = admixtures.index
    xtick_poss = bar_centers
    axes.set_xticklabels(xtick_labels, rotation=45, horizontalalignment='right')
    axes.set_xticks(xtick_poss)

    axes.set_xlim((0, bar_edges[-1]))
    axes.set_ylim((0, 1))

    if savefig:
        fig.tight_layout()
        fig.savefig(str(plot_path))

    return {'sample_order': sorted_samples}


def _plot_admixture_compositions_per_pop(admixtures, plot_path, pops_for_samples,
                                         pop_order=None):

    admixtures_per_pop = _calc_admixtures_per_pop(admixtures, pops_for_samples)
    _plot_admixtures(admixtures_per_pop, plot_path, sample_order=pop_order)


def plot_admixture_compositions_per_sample(results, out_dir,
                                           pops_for_samples=None,
                                           pop_order=None):
    for k, result in results.items():
        plot_path = out_dir / f'sample_composition.k_{k}.svg'
        _plot_admixture_compositions_per_sample(result['admixtures'], plot_path,
                                                pops_for_samples=pops_for_samples,
                                                pop_order=pop_order)


def plot_admixture_compositions_per_pop(results, out_dir, pops_for_samples,
                                        pop_order=None):
    for k, result in results.items():
        plot_path = out_dir / f'pop_composition.k_{k}.svg'
        _plot_admixture_compositions_per_pop(result['admixtures'], plot_path,
                                             pops_for_samples=pops_for_samples,
                                             pop_order=pop_order)


if __name__ == '__main__':
    prior = 'simple'

    results = parse_faststructure_results(prior)

    out_dir = config.FASTSTRUCTURE_PLOT_DIR / prior

    os.makedirs(out_dir, exist_ok=True)

    plot_path = out_dir / 'marginal_likelihoods.svg'
    plot_marginal_likelihoods(results, plot_path)

    passports = get_sample_passports()
    pops_for_samples = get_classifications_for_classification_key_path(passports, config.RANK2)

    pop_order = ['sp_pe_desert', 'sp_ec_n_wet_forest', 'slc_co', 'slc_ca', 'slc_mx',
                 'sll_mx', 'sll_vint_small',
                 'sp_pe_n_hills', 'sp_pe_n_inter-andean',
                 'sp_pe_desert_x_sp_pe_n_inter_andean',
                 'sp_pe_x_sp_ec',
                 'sp_ec_s_dry_forest',
                 'slc_ec_c_800m', 'slc_ec_n_600m', 'slc_ec_guayaquil', 'slc_ec_s_1000m',
                 'slc_pe_n_400m',
                 'slc_mx_sinaloa',
                 'slc_pe_c',
                 'slc_pe_n_1000m',
                 'slc_world', 'sll_vint',
                 'sll_moneymaker_ailsacraig', 'sll_oldbrooks', 'sll_oldgerman',
                 'sll_modern',
                 'sp_x_sl_cherry_cultivars',
                 'sp_x_sl', 
                 None
                 ]

    plot_admixture_compositions_per_pop(results, out_dir,
                                        pops_for_samples, pop_order=pop_order)

    plot_admixture_compositions_per_sample(results, out_dir,
                                           pops_for_samples, pop_order=pop_order)
