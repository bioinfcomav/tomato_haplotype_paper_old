
import config

from pprint import pprint
from collections import defaultdict
import os

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation import GT_FIELD
from variation.variations import VariationsH5, VariationsArrays

from passport import get_sample_passports
from pop_building import get_pops
import colors
from rarefaction_vars import keep_variations_variable_in_samples, NotEnoughSamplesError
from haplo import generate_sample_haplos_along_genome, create_haplo_id
from haplo_auto_classification import detected_outliers_and_classify_haplos
from rarefaction_vars import calc_pop_stats


def calc_haplo_diversity_indexes(haplos, min_num_segregating_variations, ploidy,
                                 allowed_missing_gts=0):
    polymorphism_threshold = 0.95

    haplos_array = haplos.values
    uniq_haplos, uniq_haplo_ids_that_correpond_to_all_haplos_idxs, haplo_counts = numpy.unique(haplos_array,
                                                                                               return_inverse=True,
                                                                                               return_counts=True,
                                                                                               axis=1)
    haplos_as_alleles_gts = uniq_haplo_ids_that_correpond_to_all_haplos_idxs

    res = {'num_uniq_haplos': uniq_haplos.shape[1]}
    haplo_freqs = haplo_counts / numpy.sum(haplo_counts)

    res['num_poly95'] = numpy.sum(haplo_freqs < 0.95)

    variations = VariationsArrays()
    haplos_shape = haplos.shape
    variations[GT_FIELD] = haplos_array.reshape((haplos_shape[0], haplos_shape[1], 1))

    pop_stats = calc_pop_stats(variations, allowed_missing_gts, ploidy=ploidy)
    res['pi'] = pop_stats['pi']
    res['num_variable_vars'] = pop_stats['num_variable_vars']
    res['num_poly95'] = pop_stats['num_poly95']

    variations_with_uniq_haplos_as_alleles = VariationsArrays()
    variations_with_uniq_haplos_as_alleles[GT_FIELD] = haplos_as_alleles_gts.reshape((1, haplos_as_alleles_gts.shape[0], 1))

    haplo_pop_stats = calc_pop_stats(variations_with_uniq_haplos_as_alleles, ploidy=ploidy,
                                     stats_to_calc=['unbiased_exp_het_mean'])

    res['haplo_diversity'] = haplo_pop_stats['unbiased_exp_het_mean']
    return res


def _keep_only_haplos_from_classes(haplos, haplo_classification, haplo_classes_to_keep,
                                   chrom, win_start):

        haplo_ids = [create_haplo_id(chrom, win_start, sample, haploid_idx) for (sample, haploid_idx) in haplos.columns]
        haplos_to_keep = [idx for idx, haplo_id in enumerate(haplo_ids) if haplo_classification.get(haplo_id) in haplo_classes_to_keep]
        return haplos.iloc[:, haplos_to_keep]


def do_rarefaction_for_population(variations, samples,
                                  win_params,
                                  rarefaction_range,
                                  num_wins_to_process=None,
                                  haplo_classification=None,
                                  haplo_classes_to_keep=None,
                                  allowed_missing_gts=0,
                                  ):
    ploidy = variations.ploidy
    min_num_haplotypes = rarefaction_range[0]

    if len(samples) * variations.ploidy < min_num_haplotypes:
        raise NotEnoughSamplesError()

    diversities = {}
    for haplos_info in generate_sample_haplos_along_genome(variations, win_params,
                                                           num_wins_to_process=num_wins_to_process,
                                                           samples=samples):

        haplos = haplos_info['sample_haplos']

        if haplo_classes_to_keep:    
            haplos = _keep_only_haplos_from_classes(haplos, haplo_classification, haplo_classes_to_keep,
                                                    haplos_info['chrom'], haplos_info['win_start'])

        
        num_haplos = haplos.shape[1]
        if num_haplos < min_num_haplotypes:
            continue

        num_haplos_range = range(min_num_haplotypes,
                                 min([num_haplos, rarefaction_range[1]]) + 1)


        for num_haplos in num_haplos_range:
            sampled_haplos = haplos.sample(n=num_haplos, axis=1, replace=False)
            pop_stats_for_this_num_haplos_and_this_win = calc_haplo_diversity_indexes(sampled_haplos,
                                                                                      min_num_segregating_variations=win_params['min_num_snp_for_window'],
                                                                                      allowed_missing_gts=allowed_missing_gts,
                                                                                      ploidy=ploidy)
            for diversity_index, value in pop_stats_for_this_num_haplos_and_this_win.items():
                try:
                    values_for_this_index = diversities[diversity_index]
                except KeyError:
                    values_for_this_index = defaultdict(list)
                    diversities[diversity_index] = values_for_this_index
                values_for_this_index[num_haplos].append(value)
    return diversities


def calc_rarefacted_haplo_diversities(variations, pops, rarefaction_range, win_params,
                                      num_wins_to_process=None,
                                      haplo_classification=None,
                                      allowed_missing_gts=0,
                                      haplo_classes_to_ignore=None):

    if haplo_classes_to_ignore is None:
        haplo_classes_to_ignore = []

    if haplo_classification:
        haplo_classes = [haplo_class for haplo_class in set(haplo_classification.values()) if haplo_class not in haplo_classes_to_ignore]
    else:
        haplo_classes = []
    haplo_classes.append(None)

    all_samples = {sample for samples in pops.values() for sample in samples}

    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        all_samples = all_samples.intersection(variations.samples)
    all_samples = sorted(all_samples)

    variations = keep_variations_variable_in_samples(variations, all_samples)

    diversities = {}
    for haplo_class in haplo_classes:

        haplo_classes_to_keep = None if haplo_class is None else [haplo_class]

        for pop, samples in pops.items():

            if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
                samples = sorted(set(samples).intersection(variations.samples))

            try:
                diversities[pop, haplo_class] = do_rarefaction_for_population(variations, samples,
                                                                              rarefaction_range=rarefaction_range,
                                                                              win_params=win_params,
                                                                              num_wins_to_process=num_wins_to_process,
                                                                              haplo_classification=haplo_classification,
                                                                              haplo_classes_to_keep=haplo_classes_to_keep,
                                                                              allowed_missing_gts=allowed_missing_gts)
            except NotEnoughSamplesError:
                continue    
    return diversities


def plot_rarefacted_diversities(rarefacted_diversities_per_pop, plot_path, pop_colors=None, only_these_pops=None, y_lims=None,
                                draw_percentiles=False):
    if only_these_pops is None:
        pops = sorted(rarefacted_diversities_per_pop.keys(), key=str)
    else:
        pops = only_these_pops

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    color_schema = colors.ColorSchema(pop_colors)

    if y_lims is None:
        y_lims = {}

    for pop in pops:
        x_values = []
        y_values = []
        if draw_percentiles:
            y05_values, y25_values, y75_values, y95_values = [], [], [], []
        else:
            mean_values = []

        try:
            rarefacted_values = rarefacted_diversities_per_pop[pop]
        except KeyError:
            continue

        if rarefacted_values is None:
            continue
        for num_haplos_in_rarefaction, diversities_for_pop in rarefacted_values.items():
            x_values.append(num_haplos_in_rarefaction)
            if draw_percentiles:
                per25, median, per75 = numpy.percentile(diversities_for_pop, [25, 50, 75])
                y_values.append(median)
                y25_values.append(per25)
                y75_values.append(per75)
            else:
                mean_value = numpy.mean(diversities_for_pop)
                y_values.append(mean_value)
        color = color_schema[pop]
        axes.plot(x_values, y_values, label=pop, color=color)
        if draw_percentiles:
            axes.fill_between(x_values, y25_values, y75_values, alpha=0.5, color=color)
    
    axes.legend()

    fig.tight_layout()
    fig.savefig(str(plot_path))


if __name__ == '__main__':

    haplo_classes_to_ignore = ['out_0', 'group_outlier']
    rarefaction_range = (10, 31)
    n_dims_to_keep = 3
    classification_config = {'thinning_dist_threshold': 0.00030,
                             'method': 'agglomerative',
                             'n_clusters': 3}
    classification_outlier_config = {'method': 'elliptic_envelope',
                                     'contamination': 0.015}

    outlier_configs = [{'method': 'isolation_forest', 'contamination': 0.070,
                        'thinning_dist_threshold': 0.0015}]
    classification_references = {'SL4.0ch01%610462%ts-554%1': 'sl',
                                 'SL4.0ch01%610462%ts-450%1': 'sp_peru',
                                 'SL4.0ch01%610462%bgv007339%1': 'sp_ecu'}


    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}
    cache_dir = config.CACHE_DIR

    debug = False
    if debug:
        num_wins_to_process = 2
    else:
        num_wins_to_process = None

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = get_sample_passports()

    main_pops = ['sp_pe' ,'sp_ec', 'slc_ma', 'slc_ec', 'slc_pe', 'sll_mx']
    slc = ['slc_ma', 'slc_ec', 'slc_pe', 'sll_mx']
    vintage_pops = ['sll_old_cultivars', 'sll_vint', 'slc_world']
    hybrid_pops = ['sll_modern', 'sp_x_sl', 'sp_x_sp']
    all_pops = main_pops + vintage_pops + hybrid_pops

    pops_descriptions = {config.RANK1: all_pops}
    pops = get_pops(pops_descriptions, passports)
    all_samples = {sample for samples in pops.values() for sample in samples}

    res = detected_outliers_and_classify_haplos(variations,
                                                win_params=win_params,
                                                num_wins_to_process=num_wins_to_process,
                                                samples_to_use=sorted(variations.samples),
                                                n_dims_to_keep=n_dims_to_keep,
                                                classification_config=classification_config,
                                                classification_outlier_config=classification_outlier_config,
                                                outlier_configs=outlier_configs,
                                                out_dir=None,
                                                classification_references=classification_references,
                                                cache_dir=cache_dir)
    haplo_classification = res['classification']
    
    rarefacted_haplo_diversities = calc_rarefacted_haplo_diversities(variations, pops, rarefaction_range,
                                                                     win_params=win_params,
                                                                     num_wins_to_process=num_wins_to_process,
                                                                     haplo_classification=haplo_classification,
                                                                     haplo_classes_to_ignore=haplo_classes_to_ignore
                                                                    )

    pop_colors = colors.CLASSIFICATION_RANK1_COLORS

    diversity_indexes = {diversity_index for diversities_for_pop in rarefacted_haplo_diversities.values() for diversity_index in diversities_for_pop.keys()}
    haplo_classes = {haplo_class for _, haplo_class in rarefacted_haplo_diversities.keys()}

    out_dir = config.RAREFACTION_HAPLO_DIR 
    os.makedirs(out_dir, exist_ok=True)

    for diversity_index in diversity_indexes:
        for haplo_class in haplo_classes:
            rarefacted_diversities_for_index = {}
            for (pop, this_haplo_class), diversities_for_pop in rarefacted_haplo_diversities.items():
                if this_haplo_class != haplo_class:
                    continue
                if not diversities_for_pop:
                    continue
                rarefacted_diversities_for_index[pop] = diversities_for_pop[diversity_index]

            
            for pop_group in [main_pops, vintage_pops, hybrid_pops, slc]:
                pops_str = '-'.join(map(str, pop_group))
                plot_path = out_dir / f'rarefaction_{diversity_index}.{pops_str}.{haplo_class}.svg'
                plot_rarefacted_diversities(rarefacted_diversities_for_index, plot_path=plot_path, pop_colors=pop_colors, only_these_pops=pop_group)
                plot_path = out_dir / f'rarefaction_{diversity_index}.{pops_str}.{haplo_class}.with_percentiles.svg'
                plot_rarefacted_diversities(rarefacted_diversities_for_index, plot_path=plot_path, pop_colors=pop_colors, only_these_pops=pop_group, draw_percentiles=True)
