
import config

from pprint import pprint
import random
from collections import defaultdict
from functools import partial
import os
import math

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5
from variation.variations.filters import SampleFilter, FLT_VARS

from passport import get_sample_passports
from pop_building import get_pops
from rarefaction_vars import keep_variations_variable_in_samples
from rarefaction_haplos import _keep_only_haplos_from_classes, calc_haplo_diversity_indexes
from diversities_vars import _calc_diversities, plot_diversities
from haplo import generate_sample_haplos_along_genome
from haplo_auto_classification import detected_outliers_and_classify_haplos
import colors
import labels
from diversities_haplos import calc_haplo_diversities
import matplotlib_support


def plot_diversities_by_pop_and_haplo_type(diversities_by_haplo_type, pop_haplo_type_combinations,
                                           diversity_index,
                                           plot_path, color_schema):

    bar_widths = 0.8
    haplo_separation = 0.5
    haplo_label_y_pos = 1.1
    haplo_label_font_size = 20

    #diversity_indexes = {diversity_index for diversities in diversities_by_haplo_type.values() for diversity_index in diversities.keys()}
    heights = [diversities_by_haplo_type[haplo_type].get(diversity_index, {}).get(pop, {}).get('mean', math.nan) for pop, haplo_type in pop_haplo_type_combinations]

    yerr = [diversities_by_haplo_type[haplo_type].get(diversity_index, {}).get(pop, {}).get('cim', math.nan) for pop, haplo_type in pop_haplo_type_combinations]
    haplo_types = [haplo_type for _, haplo_type in pop_haplo_type_combinations]
    colors = [color_schema[haplo_type] for haplo_type in haplo_types]

    x_values = numpy.arange(len(pop_haplo_type_combinations))

    haplo_type_idxs = []
    current_haplo_type = None
    for haplo_type in haplo_types:
        if current_haplo_type is None:
            current_haplo_type = haplo_type
            current_haplo_type_idx = 0
        if current_haplo_type is not None and current_haplo_type != haplo_type:
            current_haplo_type = haplo_type
            current_haplo_type_idx += 1
        haplo_type_idxs.append(current_haplo_type_idx)
    haplo_type_idxs = numpy.array(haplo_type_idxs)

    x_values = x_values + (haplo_type_idxs * haplo_separation)

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    axes.bar(x_values, heights, yerr=yerr, color=colors, width=bar_widths)
    x_labels = [labels.LABELS[pop] for pop, _ in pop_haplo_type_combinations]
    axes.set_xticklabels(x_labels, rotation=45, ha='right')
    axes.set_xticks(x_values)

    if diversity_index == 'mean_num_alleles':
        axes.set_ylim((1, axes.get_ylim()[1]))

    x_values_for_each_haplo_type = defaultdict(list)
    for haplo_type, x_value in zip(haplo_types, x_values):
        x_values_for_each_haplo_type[haplo_type].append(x_value)
    for haplo_type, x_values in x_values_for_each_haplo_type.items():
        x_pos = numpy.mean(x_values)
        y_pos = numpy.nanmax(heights) * haplo_label_y_pos
        axes.text(x_pos, y_pos, labels.HAPLO_LABELS[haplo_type],
                    fontsize=haplo_label_font_size,
                    color=color_schema[haplo_type],
                    horizontalalignment='center')

    matplotlib_support.set_axes_background(axes)
    axes.set_ylabel(labels.DIVERSITY_INDEX_LABELS[diversity_index])

    fig.tight_layout()
    fig.savefig(str(plot_path))


if __name__ == '__main__':

    debug = False
    
    percent_indiv_used = 0.75
    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}
    cache_dir = config.CACHE_DIR

    if debug:
        num_wins_to_process = 100
        num_repeats = 5
    else:
        num_wins_to_process = None
        num_repeats = 100

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = get_sample_passports()

    main_pops = ['sp_pe' ,'sp_ec', 'slc_ma', 'slc_ec', 'slc_pe', 'sll_mx']
    vintage_pops = ['sll_old_cultivars', 'sll_vint', 'slc_world']
    hybrid_pops = ['sll_modern', 'sp_x_sl']
    all_pops = main_pops + vintage_pops + hybrid_pops

    pops_descriptions = {config.RANK1: all_pops}
    pops = get_pops(pops_descriptions, passports)

    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        pops = {pop:sorted(set(samples).intersection(variations.samples))  for pop, samples in pops.items()}

    num_samples_per_pop ={pop: len(samples) for pop, samples in pops.items()}

    num_samples = round(min(num_samples_per_pop.values()) * percent_indiv_used)

    res = detected_outliers_and_classify_haplos(variations,
                                                win_params=win_params,
                                                num_wins_to_process=None,
                                                samples_to_use=variations.samples,
                                                n_dims_to_keep=config.N_DIMS_TO_KEEP,
                                                classification_config=config.CLASSIFICATION_CONFIG,
                                                classification_outlier_config=config.CLASSIFICATION_OUTLIER_CONFIG,
                                                outlier_configs=config.OUTLIER_CONFIGS,
                                                classification_references=config.CLASSIFICATION_REFERENCES,
                                                supervised_classification_config=config.SUPERVISED_CLASSIFICATION_CONFIG,
                                                outliers_return_aligned_pcoas=False,
                                                only_outliers = False,
                                                cache_dir=cache_dir)
    haplo_classification = res['classification']

    pop_haplo_type_combinations = [('sp_pe', 'sp_peru'), ('slc_pe', 'sp_peru'),
                                   ('sp_ec', 'sp_ecu'), ('slc_ec', 'sp_ecu'),
                                   ('slc_ma', 'sl'), ('slc_ec', 'sl'), ('slc_pe', 'sl'), ('sll_mx', 'sl'), ('sll_vint', 'sl')
                                   ]
    pop_haplo_type_combinations = [('sp_pe', 'sp_peru'), ('sp_ec', 'sp_peru'), ('slc_pe', 'sp_peru'),
                                   ('sp_pe', 'sp_ecu'), ('sp_ec', 'sp_ecu'), ('slc_ec', 'sp_ecu'),
                                   ('slc_ma', 'sl'), ('slc_ec', 'sl'), ('slc_pe', 'sl'), ('sll_mx', 'sl'), ('sll_vint', 'sl')
                                   ]


    pops_per_haplo_type = defaultdict(list)
    for (pop, haplo_type) in pop_haplo_type_combinations:
        pops_per_haplo_type[haplo_type].append(pop)

    diversities_by_haplo_type = {}
    for haplo_type, pop_names in pops_per_haplo_type.items():
        pops_for_this_haplo_type = {pop: pops[pop] for pop in pop_names}
        diversities = calc_haplo_diversities(variations, pops, num_samples, num_repeats=num_repeats,
                                             win_params=win_params,
                                             haplo_classification=haplo_classification,
                                             haplo_classes_to_keep=[haplo_type],
                                             num_wins_to_process=num_wins_to_process)
        diversities_by_haplo_type[haplo_type] =diversities

    diversity_index = 'haplo_diversity'

    out_dir = config.FIGURES_DIR
    plot_path = out_dir / 'fig5.svg'
    os.makedirs(out_dir, exist_ok=True)

    color_schema = colors.ColorSchema(colors.HAPLO_COLORS)

    plot_diversities_by_pop_and_haplo_type(diversities_by_haplo_type,
                                           pop_haplo_type_combinations,
                                           diversity_index=diversity_index,
                                           plot_path=plot_path,
                                           color_schema=color_schema)
