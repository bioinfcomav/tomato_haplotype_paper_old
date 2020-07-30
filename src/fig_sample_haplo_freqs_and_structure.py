
import config

from collections import Counter, defaultdict
import math

import numpy
import pandas

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from variation.variations import VariationsH5

import matplotlib_support
import passport
import pop_building
from haplo_auto_classification import detected_outliers_and_classify_haplos
import haplo
from faststructure_plot import plot_admixtures
from faststructure_parse_results import parse_faststructure_results
import colors
import labels


def dict_dict_to_dframe(dict_dict, fill_value=math.nan):
    index = sorted(dict_dict.keys())
    columns = sorted({key for dict_ in dict_dict.values() for key in dict_})

    array = numpy.full((len(index), len(columns)), fill_value)

    for index_idx, index_key in enumerate(index):
        for column_idx, column_key in enumerate(columns):
            array[index_idx, column_idx] = dict_dict.get(index_key, {}).get(column_key, fill_value)

    return pandas.DataFrame(array, index=index, columns=columns)


def calc_haplo_composition_per_sample(variations, samples_to_use, pops, cache_dir):

    res = detected_outliers_and_classify_haplos(variations,
                                                win_params={'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                                                            'win_size': config.HAPLO_WIN_SIZE},
                                                num_wins_to_process=None,
                                                samples_to_use=samples_to_use,
                                                n_dims_to_keep=config.N_DIMS_TO_KEEP,
                                                classification_config=config.CLASSIFICATION_CONFIG,
                                                classification_outlier_config=config.CLASSIFICATION_OUTLIER_CONFIG,
                                                outlier_configs=config.OUTLIER_CONFIGS,
                                                out_dir=config.HAPLO_PCOA_DIR,
                                                pops=pops,
                                                outliers_return_aligned_pcoas=False,
                                                only_outliers=False,
                                                classification_references=config.CLASSIFICATION_REFERENCES,
                                                supervised_classification_config=config.SUPERVISED_CLASSIFICATION_CONFIG,
                                                cache_dir=cache_dir)
    haplo_classification = res['classification']

    counts = defaultdict(Counter)
    for haplo_id, haplo_class in haplo_classification.items():
        sample = haplo.parse_haplo_id(haplo_id)[2]
        counts[sample][haplo_class] += 1

    freqs = {}
    for sample, sample_counts in counts.items():
        freqs[sample] = {haplo_class: counts_ / sum(sample_counts.values()) for haplo_class, counts_ in sample_counts.items()}

    return dict_dict_to_dframe(freqs, fill_value=0.)


if __name__ == '__main__':

    cache_dir = config.CACHE_DIR
    dist_threshold = config.CLASSIFICATION_CONFIG['thinning_dist_threshold']
    sample_for_haplo_pca = 'cervil'

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    plot_path = config.FIG_ACC_HAPLO_FREQS_AND_STRUCTURE

    sample_passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, sample_passports)
    pops_for_samples = {sample: labels.LABELS[pop] for pop, samples in pops.items() for sample in samples}
    samples_to_use = {sample for samples in pops.values() for sample in samples}

    pop_order = ['sp_pe', 'sp_montane', 'sp_x_sp', 'sp_ec',
                 'slc_ec', 'slc_co', 'slc_ma', 'slc_pe', 'slc_world',
                 'sll_mx', 'sll_vint', 'sll_modern',
                 'sp_x_sl', None]
    pop_order = [labels.LABELS[pop] for pop in pop_order]

    fig = Figure((40, 10))
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes1 = fig.add_subplot(211)

    composition_classes_colors = colors.HAPLO_COLORS.copy()
    composition_classes_colors[2] = colors.HAPLO_COLORS['sp_peru']
    composition_classes_colors[1] = colors.HAPLO_COLORS['sp_ecu']
    composition_classes_colors[0] = colors.HAPLO_COLORS['sl']

    haplo_freqs_per_sample = calc_haplo_composition_per_sample(variations, samples_to_use, pops, cache_dir)
    res = plot_admixtures(haplo_freqs_per_sample, axes=axes1, pop_order=pop_order, pops_for_samples=pops_for_samples, composition_classes_colors=composition_classes_colors)

    prior = 'simple'
    k= 3
    admixtures = parse_faststructure_results(prior)[k]['admixtures']
    axes2 = fig.add_subplot(212)
    plot_admixtures(admixtures, axes=axes2, sample_order=res['sample_order'], composition_classes_colors=composition_classes_colors)
    matplotlib_support.turn_off_x_axis(axes1)

    axes1.set_ylabel('Haplotype freq.')
    axes2.set_ylabel('Ancestral population freq.')

    fig.tight_layout()
    fig.savefig(plot_path)
