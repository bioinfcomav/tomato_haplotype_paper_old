
import config

from collections import defaultdict
import os

import pandas
import numpy

from variation.variations import VariationsH5

from passport import get_sample_passports
from pop_building import get_pops
from haplo_auto_classification import (detected_outliers_and_classify_haplos,
                                       calc_haplo_pop_composition_freq,
                                       calc_haplo_sample_composition_freq)
from haplo import get_pop_classification_for_haplos
from faststructure_parse_results import parse_faststructure_results
from plot import plot_table_classification_comparison

# TODO
# plot compositions and structure one on top of the other
# scatter plot of haplo composition sp_per with structure pop 1


def _str_or_none(value):
    if value is None:
        return 'Z' * 10
    else:
        return value


def _freqs_dict_to_dframe(composition_freq_dict):
    haplo_kinds = sorted({haplo_kind for composition in composition_freq_dict.values() for haplo_kind in composition.keys()})

    samples = sorted(composition_freq_dict.keys(), key=_str_or_none)

    freqs = defaultdict(list)
    for sample in samples:
        for haplo_kind in haplo_kinds:
            freqs[haplo_kind].append(composition_freq_dict[sample].get(haplo_kind, 0))
    dframe = pandas.DataFrame(freqs, index=map(str, samples))
    return dframe


def get_haplo_pop_composition_dframe(haplo_classification, pops):
    pop_composition_freqs = calc_haplo_pop_composition_freq(pops, haplo_classification)
    return _freqs_dict_to_dframe(pop_composition_freqs)


def get_haplo_sample_composition_dframe(haplo_classification):
    sample_composition_freqs = calc_haplo_sample_composition_freq(haplo_classification)
    return _freqs_dict_to_dframe(sample_composition_freqs)


def classify_sample_according_to_freqs(sample_composition):
    kinds = list(sample_composition.columns)

    sample_classification = [freqs.idxmax() for _, freqs in sample_composition.iterrows()]
    sample_classification = pandas.Series(sample_classification, index= sample_composition.index)
    return sample_classification


def plot_haplo_faststructure_table_comparison(haplo_classification, faststructure_results, out_dir):
    for result in faststructure_results.values():
        plot_path = out_dir / f'table_haplo_based_classification_vs_structure_classification_k_{result["k"]}.svg'
        structure_classification = classify_sample_according_to_freqs(result['admixtures'])
        plot_table_classification_comparison(haplo_classification,
                                             structure_classification,
                                             plot_path,
                                             x_label='Haplotype classification',
                                             y_label='Structure classification')


if __name__ == '__main__':

    debug = False

    num_wins_to_process = None
    cache_dir = config.CACHE_DIR
    only_outliers = False
    outliers_return_aligned_pcoas = False

    prior = 'simple'

    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = get_sample_passports()

    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, passports)

    samples_to_use = {sample for samples in pops.values() for sample in samples}
    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        samples_to_use = samples_to_use.intersection(variations.samples)
    samples_to_use = sorted(samples_to_use)

    res = detected_outliers_and_classify_haplos(variations,
                                                win_params=win_params,
                                                num_wins_to_process=num_wins_to_process,
                                                samples_to_use=samples_to_use,
                                                n_dims_to_keep=config.N_DIMS_TO_KEEP,
                                                classification_config=config.CLASSIFICATION_CONFIG,
                                                classification_outlier_config=config.CLASSIFICATION_OUTLIER_CONFIG,
                                                outlier_configs=config.OUTLIER_CONFIGS,
                                                classification_references=config.CLASSIFICATION_REFERENCES,
                                                supervised_classification_config=config.SUPERVISED_CLASSIFICATION_CONFIG,
                                                cache_dir=cache_dir)
    haplo_classification = res['classification']

    sample_haplo_composition = get_haplo_sample_composition_dframe(haplo_classification)

    haplo_sample_classification = classify_sample_according_to_freqs(sample_haplo_composition)
    faststructure_results = parse_faststructure_results(prior)

    out_dir = config.FASTSTRUCTURE_VS_HAPLO_PLOT_DIR / prior
    os.makedirs(out_dir, exist_ok=True)

    print('TODO classification threshold')

    plot_haplo_faststructure_table_comparison(haplo_sample_classification,
                                              faststructure_results, out_dir)