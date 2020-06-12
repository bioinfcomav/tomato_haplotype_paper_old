
import config

from pprint import pprint
import random
from collections import defaultdict
import os
from functools import partial
import math

import numpy

from scipy.stats import sem as calc_standard_error_of_the_mean
from scipy.stats import t as student_t

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5
from variation.variations.filters import VariableAndNotAllMissing, SampleFilter, FLT_VARS

from passport import get_sample_passports
from pop_building import get_pops
import colors
from rarefaction_vars import (keep_variations_variable_in_samples,
                              calc_pop_stats)
import colors


def _calc_diversities(variations, pops, num_samples, diversity_function,
                      num_repeats=100,
                      allowed_missing_gts=0, confidence=0.95):

    all_samples = {sample for samples in pops.values() for sample in samples}

    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        all_samples = all_samples.intersection(variations.samples)
    all_samples = sorted(all_samples)

    variations = keep_variations_variable_in_samples(variations, all_samples)

    diversities = {}
    for _ in range(num_repeats):
        for pop, samples in pops.items():
            samples_for_this_iter = random.sample(samples, num_samples)
            variations_for_this_iter = SampleFilter(samples_for_this_iter)(variations)[FLT_VARS]

            pop_stats = diversity_function(variations_for_this_iter,
                                           allowed_missing_gts=allowed_missing_gts)
            for diversity_index, value in pop_stats.items():

                if isinstance(value, (numpy.int64, numpy.float64, int, float)):
                    if diversity_index not in diversities:
                        diversities[diversity_index] = defaultdict(list)

                    diversities[diversity_index][pop].append(value)

    res = {}
    for diversity_index, diversities_for_pops in diversities.items():
        res[diversity_index] = {}
        for pop, diversity_values in diversities_for_pops.items():

            diversity_values = numpy.array(diversity_values)
            diversity_values = diversity_values[numpy.logical_not(numpy.isnan(diversity_values))]
            if not diversity_values.size:
                continue

            sem = calc_standard_error_of_the_mean(diversity_values)
            num_values = len(diversity_values)
            confidence_interval_of_the_mean = sem * student_t.ppf((1 + confidence) / 2, num_values - 1)
            res[diversity_index][pop] = {'mean': numpy.nanmean(diversity_values),
                                         'std': numpy.nanstd(diversity_values),
                                         'sem': sem,
                                         'cim': confidence_interval_of_the_mean}

    return res


def calc_var_diversities(variations, pops, num_samples, num_repeats=100,
                         allowed_missing_gts=0, confidence=0.95):


    diversity_function = partial(calc_pop_stats,
                                 allowed_missing_gts=allowed_missing_gts)
    
    return _calc_diversities(variations, pops, num_samples,
                             diversity_function=diversity_function,
                             num_repeats=num_repeats,
                             allowed_missing_gts=allowed_missing_gts, confidence=confidence)


def plot_diversities(diversities, out_dir, pop_order=None, color_schema=None):

    for diversity_index, diversities_for_index in diversities.items():

        if not diversities_for_index:
            continue

        fname = diversity_index.replace('/', '-')
        plot_path = out_dir / f'{fname}.svg'
        
        if pop_order is None:
            pop_order = sorted(diversities_for_index.keys())

        heights = [diversities_for_index.get(pop, {}).get('mean', math.nan) for pop in pop_order]
        x_values = numpy.arange(len(pop_order))

        yerr = [diversities_for_index.get(pop, {}).get('cim', math.nan) for pop in pop_order]

        if color_schema is not None:
            color = [color_schema[pop] for pop in pop_order]

        fig = Figure()
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)

        axes.bar(x_values, heights, yerr=yerr, color=color)
        axes.set_xticklabels(pop_order, rotation=45, ha='right')
        axes.set_xticks(x_values)

        if diversity_index == 'mean_num_alleles':
            axes.set_ylim((1, axes.get_ylim()[1]))

        fig.tight_layout()
        fig.savefig(str(plot_path))


if __name__ == '__main__':

    num_repeats = 100
    percent_indiv_used = 0.75

    vars_path = config.WORKING_H5
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
    pprint(num_samples_per_pop)
    num_samples = round(min(num_samples_per_pop.values()) * percent_indiv_used)

    diversities = calc_var_diversities(variations, pops, num_samples, num_repeats=num_repeats)

    out_dir = config.DIVERSITIES_VAR_DIR
    os.makedirs(out_dir, exist_ok=True)

    color_schema = colors.ColorSchema(colors.CLASSIFICATION_RANK1_COLORS)

    plot_diversities(diversities, out_dir, pop_order=all_pops, color_schema=color_schema)
