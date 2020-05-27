
import config

from pprint import pprint
from collections import defaultdict, Counter
import hashlib, pickle
import math

import numpy
import pandas
from scipy.spatial import distance

from sklearn.cluster import DBSCAN, AgglomerativeClustering, OPTICS
from sklearn.neighbors import LocalOutlierFactor
from sklearn.ensemble import IsolationForest
from sklearn.covariance import EllipticEnvelope
from sklearn.neighbors import KNeighborsClassifier

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5

from passport import get_sample_passports
from pop_building import get_pops
from haplo_pca import do_pcoas_along_the_genome, stack_aligned_pcas_projections
from procrustes import align_pcas_using_procrustes
from pca import write_curlywhirly_file
from haplo import parse_haplo_id, get_pop_classification_for_haplos
from haplo_pca_plotting import calc_ellipsoids, plot_hist2d_in_axes, plot_classifications
from util import dict_to_str
from haplo_auto_classification import detected_outliers_and_classify_haplos, thin_close_dots

HAPLO_PCOAS_X_LIMS = (-0.08, 0.03)
HAPLO_PCOAS_Y_LIMS = (-0.05, 0.03)


def get_aligned_pcoas_thinned_df(variations, dist_threshold, win_params, num_wins_to_process,
                                 samples_to_use, n_dims_to_keep,
                                 classification_config,
                                 classification_outlier_config,
                                 outlier_configs, out_dir, pops,
                                 outliers_return_aligned_pcoas,
                                 only_outliers,
                                 classification_references,
                                 supervised_classification_config,
                                 cache_dir=None):

    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += 'num_variations' + str(variations.num_variations)
        key += 'dist_threshold' + str(dist_threshold)
        key += 'win_params' + dict_to_str(win_params)
        key += 'num_wins_to_process' + str(num_wins_to_process)
        key += 'samples_to_use' + ','.join(sorted(samples_to_use))
        key += 'n_dims_to_keep' + str(n_dims_to_keep)
        key += 'outlier_config' + str([dict_to_str(outlier_config) for outlier_config in outlier_configs])
        key += 'pops' + str(pops)
        key += 'outliers_return_aligned_pcoas' + str(outliers_return_aligned_pcoas)
        key += 'only_outliers' + str(only_outliers)
        key += 'classification_references' + str(classification_references)
        key += 'supervised_classification_config' + str(supervised_classification_config)
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('thinned_aligned_haplos_pcoas_' + key + '.pickle')
        if cache_path.exists():
            return pickle.load(cache_path.open('rb'))

    res = detected_outliers_and_classify_haplos(variations,
                                                win_params=win_params,
                                                num_wins_to_process=num_wins_to_process,
                                                samples_to_use=samples_to_use,
                                                n_dims_to_keep=n_dims_to_keep,
                                                classification_config=classification_config,
                                                classification_outlier_config=classification_outlier_config,
                                                outlier_configs=outlier_configs,
                                                out_dir=out_dir,
                                                pops=pops,
                                                outliers_return_aligned_pcoas=outliers_return_aligned_pcoas,
                                                only_outliers=only_outliers,
                                                classification_references=classification_references,
                                                supervised_classification_config=supervised_classification_config,
                                                cache_dir=cache_dir)
    aligned_pcoas_df = res['aligned_pcoas_df']

    res3 = thin_close_dots(aligned_pcoas_df,
                           dist_threshold=dist_threshold)
    thinned_dots = res3['thinned_dots']
    print(thinned_dots)

    res2 = {'thinned_aligned_pcoas_df': thinned_dots, 'aligned_pcoas_df': aligned_pcoas_df,
            'classification': res['classification'], 'outlier_classes': res['outlier_classes']}
    if cache_dir:
        pickle.dump(res2, cache_path.open('wb'))

    return res2


def plot_thinned_haplos(variations, axes, dist_threshold, samples_to_use, pops, cache_dir):

    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}

    res = get_aligned_pcoas_thinned_df(variations,
                                       dist_threshold=dist_threshold,
                                       win_params=win_params,
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

    outlier_classes = res['outlier_classes']
    classification = res['classification']
    aligned_pcoas_df = res['aligned_pcoas_df']

    pop_classification = get_pop_classification_for_haplos(aligned_pcoas_df.index, pops)

    ellipsoids = calc_ellipsoids(classification, aligned_pcoas_df, classes_to_ignore=outlier_classes,
                                 scale=1.5)
    res = plot_hist2d_in_axes(res['thinned_aligned_pcoas_df'], axes,
                              x_lims=HAPLO_PCOAS_X_LIMS, y_lims=HAPLO_PCOAS_Y_LIMS,
                              ellipsoids=ellipsoids)
    return {'hist2d_result': res}


def plot_thinned_haplos2(variations, axes, dist_threshold, samples_to_use, pops,
                         cache_dir, alpha=None, haplo_colors=None):

    if haplo_colors is None:
        haplo_colors = {}

    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}

    res = get_aligned_pcoas_thinned_df(variations,
                                       dist_threshold=dist_threshold,
                                       win_params=win_params,
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

    outlier_classes = res['outlier_classes']
    classification = res['classification']
    aligned_pcoas_df = res['aligned_pcoas_df']
    aligned_thinned_pcoas_df = res['thinned_aligned_pcoas_df']

    haplo_classes = sorted(set(classification.values()))

    for haplo_class in haplo_classes:
        mask = [classification[haplo_id] == haplo_class for haplo_id in aligned_thinned_pcoas_df.index]
        x_values = aligned_thinned_pcoas_df.values[mask, 0]
        y_values = aligned_thinned_pcoas_df.values[mask, 1]
        color = haplo_colors.get(haplo_class, None)
        axes.scatter(x_values, y_values, label=haplo_class, alpha=alpha, color=color)


if __name__ == '__main__':

    cache_dir = config.CACHE_DIR
    dist_threshold = config.CLASSIFICATION_CONFIG['thinning_dist_threshold'] * 4

    out_dir = config.HAPLO_PCOA_DIR

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')
    
    passports = get_sample_passports()

    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, passports)

    samples_to_use = {sample for samples in pops.values() for sample in samples}
    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        samples_to_use = samples_to_use.intersection(variations.samples)
    samples_to_use = sorted(samples_to_use)

    plot_path = out_dir / 'thinned_pcoas_along_the_genome.hist_2d.svg'
    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    plot_thinned_haplos(variations, axes=axes, dist_threshold=dist_threshold,
                        samples_to_use=samples_to_use, pops=pops, cache_dir=cache_dir)
    fig.tight_layout()
    fig.savefig(str(plot_path))
