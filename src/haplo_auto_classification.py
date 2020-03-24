
import config

from pprint import pprint
from collections import defaultdict

import pandas

from sklearn.cluster import DBSCAN, AgglomerativeClustering, OPTICS

from variation.variations import VariationsH5

from passport import get_sample_passports
from pop_building import get_pops
from haplo_pca import do_pcoas_along_the_genome, stack_aligned_pcas_projections
from procrustes import align_pcas_using_procrustes
from pca import write_curlywhirly_file


def classify_haplo_pcoas(aligned_pcoas_df):

    method = 'agglomerative'

    if method == 'dbscan':
        clusterer = DBSCAN(eps=0.009)
    elif method == 'agglomerative':
        clusterer = AgglomerativeClustering(n_clusters=3)
    elif method == 'agglomerative2':
        clusterer = AgglomerativeClustering(n_clusters=None,
                                            distance_threshold=0.3)
    elif method == 'optics':
        clusterer = OPTICS(min_samples=0.1, max_eps=0.1)

    clustering = clusterer.fit(aligned_pcoas_df.values)

    classification = pandas.Series(clustering.labels_, index=aligned_pcoas_df.index)

    return {'classification_per_haplo_id': classification}


if __name__ == '__main__':

    debug = True

    if debug:
        num_wins_to_process = 4
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

    pops_for_samples = {sample: pop for pop, samples in pops.items() for sample in samples}
    pop_classification = {}
    for haplo_id_str in aligned_pcoas_df.index:
        sample = haplo_id_str.split('%')[2]
        pop_classification[haplo_id_str] = pops_for_samples[sample]

    res = classify_haplo_pcoas(aligned_pcoas_df)

    haplotype_clases = set(res['classification_per_haplo_id'].values)
    print(haplotype_clases)
    print('Number of haplotype classes: ', len(haplotype_clases))

    categories = {'population': pop_classification,
                  'classification': res['classification_per_haplo_id']}
    out_dir = config.HAPLO_PCOA_DIR
    out_dir.mkdir(exist_ok=True)
    path = out_dir / 'pcoas_along_the_genome.curly'
    write_curlywhirly_file(aligned_pcoas_df, path,
                           categories=categories)
