
import config

import math

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5

import passport
import pop_building
from haplo_auto_classification import detected_outliers_and_classify_haplos
from haplo_auto_classification_plot_classification_indexes import calculate_calinski_harabasz_score
from faststructure_parse_results import parse_faststructure_results
import matplotlib_support

if __name__ == '__main__':

    cache_dir = config.CACHE_DIR

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()

    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, passports)

    samples_to_use = {sample for samples in pops.values() for sample in samples}
    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        samples_to_use = samples_to_use.intersection(variations.samples)
    samples_to_use = sorted(samples_to_use)

    n_clusters_range = range(2, 9)
    classes_to_exclude = ['out_0', 'not_classified']
    res = calculate_calinski_harabasz_score(n_clusters_range, variations=variations,
                                            classes_to_exclude=classes_to_exclude,   
                                            samples_to_use=samples_to_use, cache_dir=cache_dir)

    prior = 'simple'
    results = parse_faststructure_results(prior)
    faststructure_dir = config.FASTSTRUCTURE_PLOT_DIR / prior
    faststructure_marginal_likelihoods = {k: res['marginal_likelihood'] for k, res in results.items()}

    colors = 'r', 'b'
    labels = ['Calinski Harabasz', 'Marginal likelihood']

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes1 = fig.add_subplot(111)
    axes1.plot(res['n_clusters'], res['calinski_harabasz_indexes'],
               color=colors[0])
    axes1.set_ylabel('Calinski Harabasz Index')
    matplotlib_support.set_axes_background(axes1)

    axes2 = axes1.twinx()
    faststructure_marginal_likelihoods = [faststructure_marginal_likelihoods.get(n_clusters, math.nan) for n_clusters in res['n_clusters']]
    axes2.plot(res['n_clusters'], faststructure_marginal_likelihoods,
               color=colors[1])
    axes1.set_xlabel('K')
    axes2.set_ylabel('Marginal likelihood')
    matplotlib_support.set_axes_background(axes2)
    axes2.spines['right'].set_color('grey')
    matplotlib_support.turn_off_grid(axes1)
    matplotlib_support.turn_off_grid(axes2)

    matplotlib_support.plot_legend(labels, colors, axes2, fontize=10, marker='_')

    plot_path = config.FIG_NUM_HAPLO_TYPES
    fig.tight_layout()
    fig.savefig(str(plot_path))
