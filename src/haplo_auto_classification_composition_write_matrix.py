
import config

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5

from passport import get_sample_passports
from pop_building import get_pops
from haplo_auto_classification import (detected_outliers_and_classify_haplos,
                                       calc_haplo_pop_composition_freq)
from haplo import get_pop_classification_for_haplos, parse_haplo_id




if __name__ == '__main__':

    debug = False

    outlier_configs = [{'method': 'isolation_forest', 'contamination': 0.070,
                        'thinning_dist_threshold': 0.0015}]
    n_dims_to_keep = 3
    classification_config = {'thinning_dist_threshold': 0.00025,
                             'method': 'agglomerative',
                             'n_clusters': 3}
    classification_outlier_config = {'method': 'elliptic_envelope',
                                     'contamination': 0.015}

    num_wins_to_process = None
    cache_dir = config.CACHE_DIR
    only_outliers = False
    outliers_return_aligned_pcoas = False

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
                                                n_dims_to_keep=n_dims_to_keep,
                                                classification_config=classification_config,
                                                classification_outlier_config=classification_outlier_config,
                                                outlier_configs=outlier_configs,
                                                out_dir=out_dir,
                                                cache_dir=cache_dir)

    out_path = out_dir / 'aligned_pcoas.csv.gz'
    res['aligned_pcoas_df'].to_csv(out_path, sep='\t', compression='gzip')
