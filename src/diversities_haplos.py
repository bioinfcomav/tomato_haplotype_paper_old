
import config

from pprint import pprint
import random
from collections import defaultdict
from functools import partial
import os

import numpy

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


def calc_haplo_pop_stats(variations,
                         win_params, num_wins_to_process=None,
                         haplo_classes_to_keep=None,
                         haplo_classification=None,
                         allowed_missing_gts=0):
    ploidy = variations.ploidy

    diversities = defaultdict(list)
    for haplos_info in generate_sample_haplos_along_genome(variations, win_params,
                                                           num_wins_to_process=num_wins_to_process):

        haplos = haplos_info['sample_haplos']

        if haplo_classes_to_keep:    
            haplos = _keep_only_haplos_from_classes(haplos, haplo_classification,
                                                    haplo_classes_to_keep,
                                                    haplos_info['chrom'], haplos_info['win_start'])
        pop_stats_for_this_win = calc_haplo_diversity_indexes(haplos,
                                                              min_num_segregating_variations=win_params['min_num_snp_for_window'],
                                                              allowed_missing_gts=allowed_missing_gts,
                                                              ploidy=ploidy)
        for diversity_index, value in pop_stats_for_this_win.items():
            diversities[diversity_index].append(value)

    mean_diversities = {diversity_index: numpy.nanmean(values) for diversity_index, values in diversities.items()}
    return mean_diversities


def calc_haplo_diversities(variations, pops, num_samples, win_params, num_repeats=100,
                           haplo_classification=None,
                           allowed_missing_gts=0, confidence=0.95,
                           num_wins_to_process=None):

    diversity_function = partial(calc_haplo_pop_stats,
                                 win_params=win_params,
                                 num_wins_to_process=num_wins_to_process,
                                 haplo_classification=haplo_classification,
                                 allowed_missing_gts=allowed_missing_gts)
    
    return _calc_diversities(variations, pops, num_samples,
                             diversity_function=diversity_function,
                             num_repeats=num_repeats,
                             allowed_missing_gts=allowed_missing_gts, confidence=confidence)

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
            pop_stats = calc_haplo_pop_stats(variations_for_this_iter,
                                             win_params=win_params,
                                             num_wins_to_process=num_wins_to_process,
                                             haplo_classification=haplo_classification,
                                             allowed_missing_gts=allowed_missing_gts)


if __name__ == '__main__':

    debug = False
    
    percent_indiv_used = 0.75
    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}
    cache_dir = config.CACHE_DIR

    if debug:
        num_wins_to_process = 10
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

    samples_to_use_for_haplo_classification = {sample for samples in pops.values() for sample in samples}
    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        samples_to_use_for_haplo_classification = samples_to_use_for_haplo_classification.intersection(variations.samples)
    samples_to_use_for_haplo_classification = sorted(samples_to_use_for_haplo_classification)

    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        pops = {pop:sorted(set(samples).intersection(variations.samples))  for pop, samples in pops.items()}

    num_samples_per_pop ={pop: len(samples) for pop, samples in pops.items()}
    pprint(num_samples_per_pop)
    num_samples = round(min(num_samples_per_pop.values()) * percent_indiv_used)

    res = detected_outliers_and_classify_haplos(variations,
                                                win_params=win_params,
                                                num_wins_to_process=None,
                                                samples_to_use=samples_to_use_for_haplo_classification,
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

    diversities = calc_haplo_diversities(variations, pops, num_samples, num_repeats=num_repeats,
                                         win_params=win_params,
                                         haplo_classification=haplo_classification,
                                         num_wins_to_process=num_wins_to_process)

    out_dir = config.DIVERSITIES_HAPLO_DIR
    os.makedirs(out_dir, exist_ok=True)

    color_schema = colors.ColorSchema(colors.CLASSIFICATION_RANK1_COLORS)

    plot_diversities(diversities, out_dir, pop_order=all_pops, color_schema=color_schema)
