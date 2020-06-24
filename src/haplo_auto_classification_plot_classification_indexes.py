
import config

from sklearn import metrics

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5

import passport
import pop_building
from haplo_auto_classification import detected_outliers_and_classify_haplos


def calculate_calinski_harabasz_score(n_clusters_range, variations,
                                      classes_to_exclude,
                                      samples_to_use, cache_dir):

    if not classes_to_exclude:
        classes_to_exclude = []

    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}

    classification_config = config.CLASSIFICATION_CONFIG

    n_clusterss = []
    clustering_indexes = []
    for n_clusters in n_clusters_range:
        classification_config['n_clusters'] = n_clusters
        res = detected_outliers_and_classify_haplos(variations,
                                                    win_params=win_params,
                                                    num_wins_to_process=None,
                                                    samples_to_use=samples_to_use,
                                                    n_dims_to_keep=config.N_DIMS_TO_KEEP,
                                                    classification_config=classification_config,
                                                    classification_outlier_config=config.CLASSIFICATION_OUTLIER_CONFIG,
                                                    outlier_configs=config.OUTLIER_CONFIGS,
                                                    out_dir=None,
                                                    pops=None,
                                                    outliers_return_aligned_pcoas=False,
                                                    only_outliers=False,
                                                    classification_references=config.CLASSIFICATION_REFERENCES,
                                                    supervised_classification_config=config.SUPERVISED_CLASSIFICATION_CONFIG,
                                                    cache_dir=cache_dir)
        classification = res['classification']

        if classes_to_exclude:
            classification = {haplo_id: klass for haplo_id, klass in classification.items() if klass not in classes_to_exclude}

        aligned_pcoas_df = res['aligned_pcoas_df']
        new_index = []
        labels = []
        for haplo_id in aligned_pcoas_df.index:
            try:
                label = classification[haplo_id]
            except KeyError:
                continue
            labels.append(label)
            new_index.append(haplo_id)

        aligned_pcoas_df = aligned_pcoas_df.reindex(new_index)
        clustering_index = metrics.calinski_harabasz_score(aligned_pcoas_df, labels)
        n_clusterss.append(n_clusters)
        clustering_indexes.append(clustering_index)

    return {'n_clusters': n_clusterss, 'calinski_harabasz_indexes': clustering_indexes}



if __name__ == '__main__':

    out_dir = config.HAPLO_PCOA_DIR
    out_dir.mkdir(exist_ok=True)

    cache_dir = config.CACHE_DIR

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()

    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, passports)

    samples_to_use = {sample for samples in pops.values() for sample in samples}
    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        samples_to_use = samples_to_use.intersection(variations.samples)
    samples_to_use = sorted(samples_to_use)

    out_dir = config.HAPLO_PCOA_DIR
    out_dir.mkdir(exist_ok=True)

    n_clusters_range = range(2, 9)
    classes_to_exclude = ['out_0', 'not_classified']
    res = calculate_calinski_harabasz_score(n_clusters_range, variations=variations,
                                            classes_to_exclude=classes_to_exclude,   
                                            samples_to_use=samples_to_use, cache_dir=cache_dir)

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    axes.scatter(res['n_clusters'], res['calinski_harabasz_indexes'])

    plot_path = out_dir / 'calinski_harabasz_clustering_index.svg'

    fig.savefig(str(plot_path))

