
import config

from pprint import pprint
from collections import defaultdict, Counter

import numpy
import pandas
from scipy.spatial import distance

from sklearn.cluster import DBSCAN, AgglomerativeClustering, OPTICS
from sklearn.neighbors import LocalOutlierFactor
from sklearn.ensemble import IsolationForest

from variation.variations import VariationsH5

from passport import get_sample_passports
from pop_building import get_pops
from haplo_pca import do_pcoas_along_the_genome, stack_aligned_pcas_projections
from procrustes import align_pcas_using_procrustes
from pca import write_curlywhirly_file


def classify_haplo_pcoas(aligned_pcoas_df, dist_threshold=None):

    if dist_threshold is None:
        pcoas_to_classify = aligned_pcoas_df
    else:
        thinning_result = thin_close_dots(aligned_pcoas_df, dist_threshold)
        pcoas_to_classify = thinning_result['thinned_dots']

    method = 'optics'

    if method == 'dbscan':
        clusterer = DBSCAN(eps=0.01)
    elif method == 'agglomerative':
        clusterer = AgglomerativeClustering(n_clusters=3)
    elif method == 'agglomerative2':
        clusterer = AgglomerativeClustering(n_clusters=None,
                                            distance_threshold=0.3)
    elif method == 'optics':
        clusterer = OPTICS(min_samples=0.1, max_eps=0.1)

    clustering = clusterer.fit(pcoas_to_classify.values)

    classification = clustering.labels_
    classification_index = pcoas_to_classify.index
    classification = pandas.Series(classification, index=classification_index)

    if dist_threshold is not None:
        closest_remaining_dot = thinning_result['closest_remaining_dot']
        classification = [classification.loc[closest_remaining_dot.get(haplo_idx, haplo_idx)] for haplo_idx in aligned_pcoas_df.index]

        classification = pandas.Series(classification, index=aligned_pcoas_df.index)

    return {'classification_per_haplo_id': classification}


def calc_ecuclidean_dist(dots1, dots2):
    sum_ = None
    for dim in range(dots1.shape[1]):

        this_dim_res = numpy.power(dots1[:, dim] - dots2[:, dim], 2)

        if sum_ is None:
            sum_ = this_dim_res
        else:
            sum_ += this_dim_res

    dists = numpy.sqrt(sum_)
    return dists


def thin_close_dots(dots_in_space_df, dist_threshold):
    
    thinned_dots = None
    thinned_dots_ids = []
    closest_remaining_dot = {}
    for row_id, row in dots_in_space_df.iterrows():
        #print(row_id)
        #print(row.index, row.values)
        if thinned_dots is None:
            thinned_dots = row.values
            thinned_dots = thinned_dots.reshape((1, thinned_dots.shape[0]))
            thinned_dots_ids.append(row_id)
        else:
            n_dots = n_dots = thinned_dots.shape[0]
            repeated_dot = numpy.tile(row.values, n_dots).reshape(n_dots, thinned_dots.shape[1])
            dists = calc_ecuclidean_dist(repeated_dot, thinned_dots)

            closest_dot_min_dist_idx = numpy.argmin(dists)
            min_dist = dists[closest_dot_min_dist_idx]

            if min_dist < dist_threshold:
                closest_dot_min_dist_id = thinned_dots_ids[closest_dot_min_dist_idx]
                closest_remaining_dot[row_id] = closest_dot_min_dist_id
            else:
                thinned_dots = numpy.vstack([thinned_dots , [row.values]])
                thinned_dots_ids.append(row_id)

    thinned_dots = pandas.DataFrame(thinned_dots, index=thinned_dots_ids)
    print('Num. haplos after thinning: ', thinned_dots.shape[0])
    return {'thinned_dots': thinned_dots,
            'closest_remaining_dot': closest_remaining_dot}


def detect_outliers(aligned_pcoas_df):
    method = 'isolation_forest'
    if method == 'lof':
        classifier = LocalOutlierFactor(n_neighbors=100)
    elif method == 'isolation_forest':
        classifier = IsolationForest(contamination=0.01, behaviour='deprecated')

    prediction = classifier.fit_predict(aligned_pcoas_df.values)
    is_outlier = prediction == -1

    outliers = set(aligned_pcoas_df[is_outlier].index)
    non_outliers = set(aligned_pcoas_df.index).difference(outliers)
    return {'outliers': outliers, 'non_outliers': non_outliers}


if __name__ == '__main__':
    debug = True

    thin_dist_threshold = 0.0005

    if debug:
        num_wins_to_process = 20
    else:
        num_wins_to_process = None

    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}

    vars_path = config.TIER1_PHASED_AND_IMPUTED_LOW_QUAL_09_MISSING_085
    variations = VariationsH5(str(vars_path), 'r')

    passports = get_sample_passports()

    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, passports)

    samples_to_use = {sample for samples in pops.values() for sample in samples}
    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        samples_to_use = samples_to_use.intersection(variations.samples)
    samples_to_use = sorted(samples_to_use)

    pcoas = do_pcoas_along_the_genome(variations, win_params,
                                      num_wins_to_process=num_wins_to_process,
                                      samples=samples_to_use, n_dims_to_keep=3)

    aligned_pcoas = list(align_pcas_using_procrustes(pcoas))

    aligned_pcoas_df = stack_aligned_pcas_projections(aligned_pcoas)

    print('total num. haplos: ', aligned_pcoas_df.shape[0])

    outlier_res = detect_outliers(aligned_pcoas_df)
    print('num. outliers: ', len(outlier_res['outliers']))
    aligned_pcoas_with_no_outliers = aligned_pcoas_df.loc[list(outlier_res['non_outliers'])]

    res = classify_haplo_pcoas(aligned_pcoas_with_no_outliers, thin_dist_threshold)
    non_outlier_classification = res['classification_per_haplo_id']

    outliers = outlier_res['outliers']
    classification_per_haplo_id = []
    for haplo_id in aligned_pcoas_df.index:
        klass = 'outlier' if haplo_id in outliers else non_outlier_classification.loc[haplo_id]
        classification_per_haplo_id.append(klass)
    classification_per_haplo_id = pandas.Series(classification_per_haplo_id,
                                                index=aligned_pcoas_df.index)

    print('classification counts')
    print(Counter(classification_per_haplo_id.values))

    pops_for_samples = {sample: pop for pop, samples in pops.items() for sample in samples}
    pop_classification = {}
    for haplo_id_str in aligned_pcoas_df.index:
        sample = haplo_id_str.split('%')[2]
        pop_classification[haplo_id_str] = pops_for_samples[sample]

    categories = {'population': pop_classification,
                  'classification': classification_per_haplo_id}
    out_dir = config.HAPLO_PCOA_DIR
    out_dir.mkdir(exist_ok=True)
    path = out_dir / 'pcoas_along_the_genome.curly'
    write_curlywhirly_file(aligned_pcoas_df, path,
                           categories=categories)
