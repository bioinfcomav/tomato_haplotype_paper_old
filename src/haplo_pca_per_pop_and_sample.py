
import config

from collections import defaultdict

import numpy
import pandas

from variation.variations import VariationsH5

from passport import get_sample_passports
from pop_building import get_pops
from haplo_pca_plotting import (plot_haplo_pcas, write_pcas_curly_file,
                                plot_pcas_per_pop, plot_pcas_per_sample,
                                calc_ellipsoids)
from haplo_auto_classification import (detected_outliers_and_classify_haplos,
                                       HAPLO_PCOAS_X_LIMS, HAPLO_PCOAS_Y_LIMS)


def get_sample_selection_criteria():

    rank1 = config.RANK1

    criteria = []
    samples_to_remove = []
    samples_to_keep = []

    return {'criteria': criteria, 'samples_to_remove': samples_to_remove,
            'samples_to_keep': samples_to_keep}


def get_haplotypes_to_exclude(path):
    haplos_to_exclude = defaultdict(list)

    if not path.exists():
        return {}

    for line in path.open('rt'):
        if line.startswith('Label'):
            continue
        line = line.strip()
        if not line:
            continue
        items = line.strip().split()[0].split('%')
        chrom = f'SL2.50ch{items[0]}'
        win_start = int(items[1])
        key = chrom, win_start
        sample = items[2]
        haploid_idx = int(items[3])
        haplos_to_exclude[key].append((sample, haploid_idx))
    return haplos_to_exclude


def get_haplotypes_to_include(path):
    return get_haplotypes_to_exclude(path)


if __name__ == '__main__':

    debug = False

    if debug:
        num_wins_to_process = 2
        cache_dir = None
    else:
        num_wins_to_process = None
        cache_dir = config.CACHE_DIR

    outlier_configs = [{'method': 'isolation_forest', 'contamination': 0.070,
                        'thinning_dist_threshold': 0.0015}]
    n_dims_to_keep = 3
    classification_config = {'thinning_dist_threshold': 0.00025,
                             'method': 'agglomerative',
                             'n_clusters': 3}
    classification_outlier_config = {'method': 'elliptic_envelope',
                                     'contamination': 0.015}
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

    out_dir = config.HAPLO_PCOA_DIR
    out_dir.mkdir(exist_ok=True)

    res = detected_outliers_and_classify_haplos(variations,
                                                win_params=win_params,
                                                num_wins_to_process=num_wins_to_process,
                                                samples_to_use=samples_to_use,
                                                n_dims_to_keep=config.N_DIMS_TO_KEEP,
                                                classification_config=config.CLASSIFICATION_CONFIG,
                                                classification_outlier_config=config.CLASSIFICATION_OUTLIER_CONFIG,
                                                outlier_configs=config.OUTLIER_CONFIGS,
                                                out_dir=out_dir,
                                                classification_references=config.CLASSIFICATION_REFERENCES,
                                                supervised_classification_config=config.SUPERVISED_CLASSIFICATION_CONFIG,
                                                cache_dir=cache_dir)
    haplo_classification = res['classification']
    aligned_pcoas_df = res['aligned_pcoas_df']

    ellipsoids = calc_ellipsoids(haplo_classification, aligned_pcoas_df, classes_to_ignore=['out_0', 'group_outlier'],
                                 scale=1.5)


    per_pop_out_dir = out_dir / 'per_pop'
    per_pop_out_dir.mkdir(exist_ok=True)
    plot_pcas_per_pop(aligned_pcoas_df, per_pop_out_dir,
                      populations=pops, ellipsoids=ellipsoids,
                      x_lims=HAPLO_PCOAS_X_LIMS, y_lims=HAPLO_PCOAS_Y_LIMS)

    per_sample_out_dir = out_dir / 'per_sample'
    per_sample_out_dir.mkdir(exist_ok=True)
    plot_pcas_per_sample(aligned_pcoas_df, per_sample_out_dir,
                         x_lims=HAPLO_PCOAS_X_LIMS, y_lims=HAPLO_PCOAS_Y_LIMS,
                         populations=pops,
                         ellipsoids=ellipsoids)
