
import config

from pprint import pprint
import itertools
from array import array
import operator

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5

from haplo import generate_dists_between_sample_haplos_along_genome
import passport
import pop_building
import labels
from haplo_auto_classification import detected_outliers_and_classify_haplos


def classify_dists(dists, haplo_classification):
    for dist in dists:
        haplo_class1 = haplo_classification[dist['haplo1_id']]
        haplo_class2 = haplo_classification[dist['haplo2_id']]
        haplo_classes = sorted([haplo_class1, haplo_class2])
        within_class = haplo_class1 == haplo_class2
        yield {'haplo_class1': haplo_class1,
               'haplo_class2': haplo_class2,
               'within_class': within_class,
               'dist': dist['dist']}


def plot_dist_hist_by_haplo_class(dists, axes, dist_range=None, num_bins=15, standarize=False):
    dists_by_class = {True: array('f'), False: array('f')}
    for dist in dists:
        dists_by_class[dist['within_class']].append(dist['dist'])

    dists_by_class = {klass: numpy.array(dists) for klass, dists in dists_by_class.items()}

    if dist_range is None:
        dist_range = (min([numpy.nanmin(dists) for dists in dists_by_class.values() if dists.size]),
                      max([numpy.nanmax(dists) for dists in dists_by_class.values() if dists.size]))

    bin_edges = numpy.linspace(dist_range[0], dist_range[1], num=num_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    for dist_class, dists in dists_by_class.items():
        if not dists.size:
            continue
        label = 'within haplotype class' if dist_class else 'between haplotype classes'
        counts, _ = numpy.histogram(dists, bins=bin_edges)
        if standarize:
            max_count = numpy.nanmax(counts)
            counts = counts / max_count

        axes.plot(bin_centers, counts, label=label)


def plot_dist_hist(dists, axes, dist_range=None, num_bins=10, standarize=False, label=None):

    dists = array('f', map(operator.itemgetter('dist'), dists))

    if not dists:
        return

    dists = numpy.array(dists)

    if dist_range is None:
        dist_range = (numpy.nanmin(dists), numpy.nanmax(dists))

    bin_edges = numpy.linspace(dist_range[0], dist_range[1], num=num_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    counts, _ = numpy.histogram(dists, bins=bin_edges)
    if standarize:
        max_count = numpy.nanmax(counts)
        counts = counts / max_count

    axes.plot(bin_centers, counts, label=label)



if __name__ == '__main__':
    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}

    vars_path = config.WORKING_H5
    variations = VariationsH5(str(vars_path), 'r')
    samples = None
    num_wins_to_process = None
    min_num_snp_for_dist = 0
    num_dists_to_process = 1000000
    dist_range = (0, 5e-5)
    standarize = True
    remove_zeros = False

    sample_passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, sample_passports)

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    imputed_variations = VariationsH5(str(vars_path), 'r')

    samples_to_use = {sample for samples in pops.values() for sample in samples}
    res = detected_outliers_and_classify_haplos(variations,
                                                win_params={'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                                                            'win_size': config.HAPLO_WIN_SIZE},
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
                                                cache_dir=config.CACHE_DIR)
    haplo_classification = res['classification']

    required_pops = ['sp_pe', 'sp_ec', 'slc_ma', 'slc_ec', 'slc_pe', 'sp_x_sl']

    out_dir = config.TMP_DIR
    out_dir.mkdir(exist_ok=True)

    if remove_zeros:
        plot_path = out_dir / f'dist_hist_no_zeros.svg'
    else:
        plot_path = out_dir / f'dist_hist.svg'

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    for pop in required_pops:
        samples = pops[pop]
        dists = generate_dists_between_sample_haplos_along_genome(variations, win_params,
                                                                min_num_snp_for_dist=min_num_snp_for_dist,
                                                                num_wins_to_process=num_wins_to_process,
                                                                samples=samples)
        #classified_dists = classify_dists(dists, haplo_classification)

        if remove_zeros:
            dists = filter(lambda x: x['dist'] > 0, dists)

        if num_dists_to_process:
            dists = itertools.islice(dists, num_dists_to_process)

        #plot_dist_hist_by_haplo_class(classified_dists, axes, dist_range=dist_range, standarize=standarize)
        plot_dist_hist(dists, axes, dist_range=dist_range, standarize=standarize, label=pop)
    axes.legend()
    fig.savefig(str(plot_path))
    