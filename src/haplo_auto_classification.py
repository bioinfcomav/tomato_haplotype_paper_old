
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

from variation.variations import VariationsH5

from passport import get_sample_passports
from pop_building import get_pops
from haplo_pca import do_pcoas_along_the_genome, stack_aligned_pcas_projections
from procrustes import align_pcas_using_procrustes
from pca import write_curlywhirly_file
from haplo import parse_haplo_id, get_pop_classification_for_haplos
from haplo_pca_plotting import calc_ellipsoids, plot_hist2d, plot_classifications
from util import dict_to_str

HAPLO_PCOAS_X_LIMS = (-0.08, 0.03)
HAPLO_PCOAS_Y_LIMS = (-0.05, 0.03)


def _classify_haplo_pcoas(aligned_pcoas_df, classification_config):

    params = classification_config.copy()

    thinning_dist_threshold = params['thinning_dist_threshold']
    del params['thinning_dist_threshold']

    if thinning_dist_threshold is None:
        pcoas_to_classify = aligned_pcoas_df
    else:
        thinning_result = thin_close_dots(aligned_pcoas_df, thinning_dist_threshold)
        pcoas_to_classify = thinning_result['thinned_dots']

    method = params['method']
    del params['method']

    if method == 'dbscan':
        clusterer = DBSCAN(**params)
    elif method == 'agglomerative':
        clusterer = AgglomerativeClustering(**params)
    elif method == 'optics':
        clusterer = OPTICS(**params)

    clustering = clusterer.fit(pcoas_to_classify.values)

    classification = clustering.labels_
    classification_index = pcoas_to_classify.index
    classification = pandas.Series(classification, index=classification_index)

    if thinning_dist_threshold is not None:
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


def _detect_outliers(aligned_pcoas_df, outlier_config):

    params = outlier_config.copy()
    if 'thinning_dist_threshold' in params:
        thinning_dist_threshold = params['thinning_dist_threshold']
        del params['thinning_dist_threshold']
    else:
        thinning_dist_threshold = None
    method = outlier_config['method']

    del params['method']

    if thinning_dist_threshold is None:
        thinned_aligned_pcoas_df = aligned_pcoas_df
    else:
        thinning_result = thin_close_dots(aligned_pcoas_df, thinning_dist_threshold)
        thinned_aligned_pcoas_df = thinning_result['thinned_dots']

    if method == 'lof':
        classifier = LocalOutlierFactor(**params)
    elif method == 'isolation_forest':
        classifier = IsolationForest(**params)
    elif method == 'elliptic_envelope':
        classifier = EllipticEnvelope(**params)
    else:
        raise ValueError(f'unkown method: {method}')

    prediction = classifier.fit_predict(thinned_aligned_pcoas_df.values)
    is_outlier = prediction == -1
    is_outlier = pandas.Series(is_outlier, index=thinned_aligned_pcoas_df.index)

    if thinning_dist_threshold is not None:
        closest_remaining_dot = thinning_result['closest_remaining_dot']
        is_outlier = [is_outlier.loc[closest_remaining_dot.get(haplo_idx, haplo_idx)] for haplo_idx in aligned_pcoas_df.index]
        is_outlier = pandas.Series(is_outlier, index=aligned_pcoas_df.index)

    outliers = set(aligned_pcoas_df[is_outlier.values].index)
    non_outliers = set(aligned_pcoas_df.index).difference(outliers)
    
    return {'outliers': outliers, 'non_outliers': non_outliers,
            'is_outlier': is_outlier}


def _detect_outlier_haplos(variations, win_params, num_wins_to_process,
                           samples_to_use, n_dims_to_keep,
                           outlier_config,
                           haplotypes_to_exclude=None,
                           aligned_pcoas_df=None,
                           cache_dir=None):
    if aligned_pcoas_df is not None and haplotypes_to_exclude is not None:
        raise ValueError('if a pcoa is provided you can not ask for haplos to be excluded in the pcoa')

    if aligned_pcoas_df is None:

        pcoas = do_pcoas_along_the_genome(variations, win_params,
                                          num_wins_to_process=num_wins_to_process,
                                          samples=samples_to_use,
                                          n_dims_to_keep=n_dims_to_keep,
                                          haplotypes_to_exclude=haplotypes_to_exclude,
                                          cache_dir=cache_dir)

        aligned_pcoas = list(align_pcas_using_procrustes(pcoas))

        aligned_pcoas_df = stack_aligned_pcas_projections(aligned_pcoas)

    print('total num. haplos: ', aligned_pcoas_df.shape[0])

    outlier_res = _detect_outliers(aligned_pcoas_df,
                                   outlier_config=outlier_config)
    print('num. outliers: ', len(outlier_res['outliers']))
    outlier_res['aligned_pcoas_df'] = aligned_pcoas_df

    return outlier_res


def collect_outlier_haplos_by_win(outlier_haplos):

    outlier_haplos_by_win = defaultdict(list)
    for haplo in set(outlier_haplos):
        chrom, win_start, sample, haploid_idx = parse_haplo_id(haplo)
        outlier_haplos_by_win[(chrom, win_start)].append((sample, haploid_idx))
    return outlier_haplos_by_win


def detect_outlier_haplos(variations, win_params, num_wins_to_process,
                          samples_to_use, n_dims_to_keep,
                          outlier_configs, out_dir, pops,
                          return_aligned_pcoas,
                          cache_dir):

    if cache_dir and return_aligned_pcoas:
        raise ValueError('cache cannot be used with aligned PCoAs')

    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += 'num_variations' + str(variations.num_variations)
        key += 'win_params' + dict_to_str(win_params)
        key += 'num_wins_to_process' + str(num_wins_to_process)
        key += 'samples_to_use' + ','.join(sorted(samples_to_use))
        key += 'n_dims_to_keep' + str(n_dims_to_keep)
        key += 'outlier_config' + str([dict_to_str(outlier_config) for outlier_config in outlier_configs])
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('outlier_haplos' + key + '.pickle')
        if cache_path.exists():
            return pickle.load(cache_path.open('rb'))

    first_aligned_pcoas_df = None
    aligned_pcoas_df = None
    outlier_classes = {}
    for idx, outlier_config in enumerate(outlier_configs):
        if aligned_pcoas_df is None:
            outlier_res = _detect_outlier_haplos(variations, win_params, num_wins_to_process,
                                                 samples_to_use, n_dims_to_keep,
                                                 outlier_config,
                                                 haplotypes_to_exclude=collect_outlier_haplos_by_win(outlier_classes.keys()),
                                                 cache_dir=cache_dir)
            aligned_pcoas_df = outlier_res['aligned_pcoas_df']
        else:
            current_outliers = set(outlier_classes.keys())
            mask = [haplo_id not in current_outliers for haplo_id in aligned_pcoas_df.index]
            aligned_pcoas_df = aligned_pcoas_df[mask]
            outlier_res = _detect_outlier_haplos(variations, win_params, num_wins_to_process,
                                                 samples_to_use, n_dims_to_keep,
                                                 outlier_config,
                                                 aligned_pcoas_df=aligned_pcoas_df)
        if first_aligned_pcoas_df is None:
            first_aligned_pcoas_df = outlier_res['aligned_pcoas_df']
        for sample, is_outlier in outlier_res['is_outlier'].iteritems():
            if is_outlier:
                outlier_classes[sample] = f'out_{idx}'

    categories = {'outliers': outlier_classes}
    if pops:
        pops_for_samples = {sample: pop for pop, samples in pops.items() for sample in samples_to_use}
        pop_classification = {}
        for haplo_id_str in first_aligned_pcoas_df.index:
            sample = haplo_id_str.split('%')[2]
            pop_classification[haplo_id_str] = pops_for_samples[sample]
        categories['population'] = pop_classification

    if out_dir is not None:
        path = out_dir / 'pcoas_along_the_genome.with_outliers.curly'
        write_curlywhirly_file(first_aligned_pcoas_df, path,
                            categories=categories)

    res = {'outlier_classes': outlier_classes}

    if return_aligned_pcoas:
        res['aligned_pcoas_df'] = first_aligned_pcoas_df

    if cache_dir:
        pickle.dump(res, cache_path.open('wb'))

    return res


def remove_outliers_from_classified_clusters(classification, aligned_pcoas_df,
                                             outlier_config):

    haplo_ids_by_class = defaultdict(set)
    for haplo_id, klass in classification.items():
        haplo_ids_by_class[klass].add(haplo_id)

    outlier_names = {klass: f'outlier_from_group_{klass}' for klass in haplo_ids_by_class.keys()}

    for klass, haplo_ids in haplo_ids_by_class.items():
        aligned_pcoas_df_for_this_class = aligned_pcoas_df.loc[haplo_ids, :]

        res = _detect_outliers(aligned_pcoas_df_for_this_class, outlier_config)
        for haplo_id in res['outliers']:
            classification[haplo_id] = outlier_names[klass]
    return {'outlier_classes': outlier_names}


def _do_supervissed_classification(aligned_pcoas_df, labelled_haplo_ids, config):

    known_haplo_ids, known_labels = list(zip(*labelled_haplo_ids.items()))

    aligned_pcoas_df_for_known = aligned_pcoas_df.reindex(known_haplo_ids)

    classifier_config = config.copy()
    del classifier_config['prob_threshold']
    del classifier_config['classifier']

    if config['classifier'] == 'kneighbors':
        classifier = KNeighborsClassifier(**classifier_config)

    classifier.fit(aligned_pcoas_df_for_known.values, known_labels)
    probs = classifier.predict_proba(aligned_pcoas_df)

    num_classes =  probs.shape[1]
    haplo_ids = numpy.array(aligned_pcoas_df.index)
    classification = {}
    for klass in range(num_classes):
        mask = probs[:, klass] >= config['prob_threshold']
        for haplo_id in haplo_ids[mask]:
            classification[haplo_id] = klass

    for haplo_id in haplo_ids:
        if haplo_id not in classification:
            classification[haplo_id] = 'not_classified'
    return {'classification': classification,
            'outlier_classes': ['not_classified']}


def classify_haplos(variations, win_params, num_wins_to_process,
                    samples, n_dims_to_keep,
                    outlier_classes,
                    classification_config,
                    supervised_classification_config,
                    outlier_config, cache_dir=None):

    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += 'num_variations' + str(variations.num_variations)
        key += 'win_params' + dict_to_str(win_params)
        key += 'num_wins_to_process' + str(num_wins_to_process)
        key += 'samples_to_use' + ','.join(sorted(samples))
        key += 'outlier_clases' + str(sorted(outlier_classes.values()))
        key += 'n_dims_to_keep' + str(n_dims_to_keep)
        key += 'clasfication_config' + dict_to_str(classification_config)
        key += 'outlier_config' + dict_to_str(outlier_config)
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('haplo_clasfication' + key + '.pickle')
        if cache_path.exists():
            return pickle.load(cache_path.open('rb'))

    outlier_haplos_by_win = collect_outlier_haplos_by_win(outlier_classes.keys())

    pcoas = do_pcoas_along_the_genome(variations, win_params,
                                      num_wins_to_process=num_wins_to_process,
                                      samples=samples, n_dims_to_keep=n_dims_to_keep,
                                      haplotypes_to_exclude=outlier_haplos_by_win,
                                      cache_dir=cache_dir)

    aligned_pcoas = list(align_pcas_using_procrustes(pcoas))

    aligned_pcoas_df = stack_aligned_pcas_projections(aligned_pcoas)

    print('total num. non outlier haplos: ', aligned_pcoas_df.shape[0])

    res = _classify_haplo_pcoas(aligned_pcoas_df,
                                classification_config)
    classification = dict(res['classification_per_haplo_id'].iteritems())
    haplo_classes = set(classification.values())

    res = remove_outliers_from_classified_clusters(outlier_config=outlier_config, 
                                                   classification=classification,
                                                   aligned_pcoas_df=aligned_pcoas_df)
    #outlier_classes = res['outlier_classes']

    classification.update(outlier_classes)

    # now we do a suppervised classification
    core_classification = {haplo_id: klass for haplo_id, klass in classification.items() if klass in haplo_classes}

    res = _do_supervissed_classification(aligned_pcoas_df, core_classification, supervised_classification_config)

    classification.update(res['classification'])
    outlier_classes = res['outlier_classes']

    res = {'classification': classification, 'aligned_pcoas_df': aligned_pcoas_df,
           'outlier_classes': outlier_classes}

    if cache_dir:
        pickle.dump(res, cache_path.open('wb'))

    return res


def rename_classification(classification, classification_references):
    mapping = {classification[haplo_id]: nice_klass_name for haplo_id, nice_klass_name in classification_references.items()}

    return {haplo_id: mapping.get(klass, klass) for haplo_id, klass in classification.items()}


def detected_outliers_and_classify_haplos(variations, win_params,
                                          num_wins_to_process,
                                          samples_to_use,
                                          n_dims_to_keep,
                                          classification_config,
                                          classification_outlier_config,
                                          outlier_configs,
                                          supervised_classification_config,
                                          out_dir=None,
                                          pops=None,
                                          outliers_return_aligned_pcoas=False,
                                          only_outliers=False,
                                          classification_references=None,
                                          cache_dir=None):
    outlier_classes = []
    # first iteration, with outliers,
    res = detect_outlier_haplos(variations, win_params=win_params,
                                num_wins_to_process=num_wins_to_process,
                                samples_to_use=samples_to_use,
                                n_dims_to_keep=n_dims_to_keep,
                                outlier_configs=outlier_configs,
                                pops=pops,
                                out_dir=out_dir,
                                cache_dir=cache_dir,
                                return_aligned_pcoas=outliers_return_aligned_pcoas)
    outlier_classification = res['outlier_classes']
    outlier_classes.extend(set(outlier_classification.values()))
    if outliers_return_aligned_pcoas:
        aligned_pcoas_df = res['aligned_pcoas_df']

    # second iteration, with no outliers
    if not only_outliers:

        res = classify_haplos(variations, win_params=win_params,
                              num_wins_to_process=num_wins_to_process,
                              samples=samples_to_use,
                              n_dims_to_keep=n_dims_to_keep,
                              outlier_classes=outlier_classification,
                              classification_config=classification_config,
                              outlier_config=classification_outlier_config,
                              supervised_classification_config=supervised_classification_config,
                              cache_dir=cache_dir)
        classification = res['classification']
        aligned_pcoas_df = res['aligned_pcoas_df']
        outlier_classes.extend(res['outlier_classes'])
    else:
        classification = outlier_classification

    if classification_references:
        classification = rename_classification(classification, classification_references)

    return {'classification': classification, 'aligned_pcoas_df': aligned_pcoas_df,
            'outlier_classes': outlier_classes}


def calc_haplo_sample_composition(haplo_classification):
    sample_haplo_composition = defaultdict(Counter)
    for haplo_id, klass in haplo_classification.items():
        sample = parse_haplo_id(haplo_id)[2]
        sample_haplo_composition[sample][klass] += 1
    return sample_haplo_composition


def calc_haplo_sample_composition_freq(haplo_classification):
    sample_haplo_composition_freqs = {}
    for sample, composition in calc_haplo_sample_composition(haplo_classification).items():
        tot_count = sum(composition.values())
        sample_haplo_composition_freqs[sample] = {klass: count / tot_count for klass, count in composition.items()}
    return sample_haplo_composition_freqs


def calc_haplo_pop_composition_freq(pops, haplo_classification):
    sample_haplo_composition_freqs = calc_haplo_sample_composition_freq(haplo_classification)

    haplo_classes = {klass for composition in sample_haplo_composition_freqs.values() for klass in composition.keys()}

    pop_composition_freqs = {}
    for pop, samples in pops.items():
        pop_compositions = defaultdict(list)
        for sample in samples:
            try:
                sample_haplo_composition = sample_haplo_composition_freqs[sample]
            except KeyError:
                if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
                    continue
                else:
                    raise

            for klass in haplo_classes:
                pop_compositions[klass].append(sample_haplo_composition.get(klass, 0))
        pop_composition_freqs[pop] = {klass: numpy.mean(freqs) for klass, freqs in pop_compositions.items()}
        assert math.isclose(sum(pop_composition_freqs[pop].values()), 1)
    return pop_composition_freqs


def calc_haplo_pop_composition_freq_dframe(pops, haplo_classification):
    pop_composition_freqs = calc_haplo_pop_composition_freq(pops, haplo_classification)
    pops = list(pop_composition_freqs.keys())
    haplo_classes = sorted({haplo_class for freqs in pop_composition_freqs.values() for haplo_class in freqs.keys()})
    freqs = numpy.empty((len(pops), len(haplo_classes)))
    for pop_idx, pop in enumerate(pops):
        for class_idx, klass in enumerate(haplo_classes):
            freqs[pop_idx, class_idx] = pop_composition_freqs[pop][klass]
    return pandas.DataFrame(freqs, index=pops, columns=haplo_classes)


def group_haploids_by_class(classification):
    classification_by_class = defaultdict(list)
    for haplo_id, klass in classification.items():
        classification_by_class[klass].append(haplo_id)
    return classification_by_class


if __name__ == '__main__':

    debug = False

    if False:
        num_wins_to_process = 100
        #num_wins_to_process = 2
        only_outliers = True
        cache_dir = None
        outliers_return_aligned_pcoas = True
    elif debug:
        num_wins_to_process = 100
        #num_wins_to_process = 20
        #num_wins_to_process = 2
        only_outliers = False
        cache_dir = config.CACHE_DIR
        cache_dir = None
        outliers_return_aligned_pcoas = False
    else:
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

    classification_config = config.CLASSIFICATION_CONFIG
    classification_config['n_clusters'] = 8

    res = detected_outliers_and_classify_haplos(variations,
                                                win_params=win_params,
                                                num_wins_to_process=num_wins_to_process,
                                                samples_to_use=samples_to_use,
                                                n_dims_to_keep=config.N_DIMS_TO_KEEP,
                                                classification_config=classification_config,
                                                classification_outlier_config=config.CLASSIFICATION_OUTLIER_CONFIG,
                                                outlier_configs=config.OUTLIER_CONFIGS,
                                                out_dir=out_dir,
                                                pops=pops,
                                                outliers_return_aligned_pcoas=outliers_return_aligned_pcoas,
                                                only_outliers=only_outliers,
                                                classification_references=config.CLASSIFICATION_REFERENCES,
                                                supervised_classification_config=config.SUPERVISED_CLASSIFICATION_CONFIG,
                                                cache_dir=cache_dir)
    classification = res['classification']
    print(len(classification))

    aligned_pcoas_df = res['aligned_pcoas_df']
    outlier_classes = res['outlier_classes']
    print('classification counts')
    pprint(Counter(classification.values()))

    pop_classification = get_pop_classification_for_haplos(aligned_pcoas_df.index, pops)

    ellipsoids = calc_ellipsoids(classification, aligned_pcoas_df, classes_to_ignore=outlier_classes,
                                 scale=1.5)

    path = out_dir / 'pcoas_along_the_genome.hist_2d.svg'
    plot_hist2d(aligned_pcoas_df, path,
                x_lims=HAPLO_PCOAS_X_LIMS, y_lims=HAPLO_PCOAS_Y_LIMS,
                ellipsoids=ellipsoids)
    path = out_dir / 'pcoas_along_the_genome.hist_2d.classifications.svg'
    haplo_classification_by_class = group_haploids_by_class(classification)
    plot_classifications(aligned_pcoas_df, haplo_classification_by_class, plot_path=path,
                         outlier_classes=outlier_classes,
                         classes_to_ignore=['out_0'],
                         x_lims=HAPLO_PCOAS_X_LIMS, y_lims=HAPLO_PCOAS_Y_LIMS,
                         ellipsoids=ellipsoids)

    categories = {'population': pop_classification,
                  'classification': classification}
    path = out_dir / 'pcoas_along_the_genome.without_outliers.curly'
    write_curlywhirly_file(aligned_pcoas_df, path,
                           categories=categories)
