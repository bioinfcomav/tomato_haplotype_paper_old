
import config

import os
import random
from array import array

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation import GT_FIELD, CHROM_FIELD, POS_FIELD
from variation.variations import VariationsH5
from variation.variations.filters import SampleFilter, FLT_VARS

import passport
import pop_building
import colors
from abba import _filter_vars
from rarefaction_haplos import (_keep_only_haplos_from_classes,
                                calc_haplo_diversity_indexes)
from haplo import generate_sample_haplos_along_genome
import matplotlib_support
import labels

'''
Caicedo

Genomic Evidence for Complex Domestication History of the Cultivated Tomato in Latin America 

We used R to calculate genome-wide reduction in average π in the descendant
populations compared with their probable ancestors, calculated as a
ratio (πancestor/πdescendant) within 10-kb windows. To avoid loss of information
for the windows without genetic diversity (π = 0) in the descendant population, 
we changed the π values for those windows to the minimum window π in the descendant
population. To detect genomic regions with high genetic differentiation,
we calculated Weir and Cockerham’s FST (Weir and Cockerham 1984), using VCFtools,
and deviation of SFS from neutrality, using SweeD (Pavlidis et al. 2013), for each 
10-kb window. For FST calculations, we excluded SNPs with MAF <2.5% to adjust for
the effect of rare alleles on FST values.'''


def calc_pi_along_the_genome_for_pop(variations, samples,
                                     win_params, allowed_missing_gts=0,
                                     num_wins_to_process=None,
                                     haplo_classification=None,
                                     haplo_classes_to_keep=None):
    ploidy = variations.ploidy

    pis = array('f')

    for haplos_info in generate_sample_haplos_along_genome(variations, win_params,
                                                        num_wins_to_process=num_wins_to_process):

        haplos = haplos_info['sample_haplos']

        if haplo_classes_to_keep:    
            haplos = _keep_only_haplos_from_classes(haplos, haplo_classification,
                                                    haplo_classes_to_keep,
                                                    haplos_info['chrom'], haplos_info['win_start'])
        if not haplos.size:
            continue

        pop_stats_for_this_win = calc_haplo_diversity_indexes(haplos,
                                                                min_num_segregating_variations=win_params['min_num_snp_for_window'],
                                                                allowed_missing_gts=allowed_missing_gts,
                                                                ploidy=ploidy)
        pis.append(pop_stats_for_this_win['pi'])

    return {'pis': numpy.array(pis)}


def calculate_pi_along_the_genome(variations, pops, num_samples,
                                  win_params, allowed_missing_gts=0,
                                  num_wins_to_process=None):

    pis_per_pop = {}
    for pop, samples in pops.items():
        samples_for_this_iter = random.sample(samples, num_samples)
        variations_for_this_pop = SampleFilter(samples_for_this_iter)(variations)[FLT_VARS]

        pis = calc_pi_along_the_genome_for_pop(variations_for_this_pop, samples_for_this_iter,
                                       win_params=win_params,
                                       allowed_missing_gts=allowed_missing_gts,
                                       num_wins_to_process=num_wins_to_process)
        pis_per_pop[pop] = pis['pis']
    return pis_per_pop


def plot_distributions_per_pop(values_per_pop, axes, range_=None, n_bins=10,
                               color_schema=None):

    if color_schema is None:
        color_schema = colors.ColorSchema()

    if range_ is None:
        min_ = min([numpy.nanmin(values) for values in values_per_pop.values()])
        max_ = max([numpy.nanmax(values) for values in values_per_pop.values()])
        range_ = (min_, max_)

    bin_edges = numpy.linspace(min_, max_, n_bins + 1)
    x_poss = (bin_edges[:-1] + bin_edges[1:]) / 2

    for pop, values in values_per_pop.items():
        counts, _ = numpy.histogram(values, bins=bin_edges)
        axes.plot(x_poss, counts, label=pop, color=color_schema[pop])


if __name__ == '__main__':

    debug = False

    percent_indiv_used = 0.75
    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}
    cache_dir = config.CACHE_DIR

    if debug:
        num_wins_to_process = 20
        num_repeats = 5
    else:
        num_wins_to_process = None
        num_repeats = 100

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()

    pop_names = ['sp_pe', 'sp_ec', 'slc_ma', 'slc_ec', 'slc_pe_n', 'slc_pe_s', 'sll_mx']

    pops_descriptions = {config.RANK1: pop_names}
    pops = pop_building.get_pops(pops_descriptions, passports)
    all_samples = [sample for samples in pops.values() for sample in samples]

    num_samples_per_pop ={pop: len(samples) for pop, samples in pops.items()}

    num_samples = round(min(num_samples_per_pop.values()) * percent_indiv_used)

    fields_to_keeep = [GT_FIELD, CHROM_FIELD, POS_FIELD]
    variations = _filter_vars(variations, all_samples, fields_to_keeep=fields_to_keeep)

    pis = calculate_pi_along_the_genome(variations, pops, num_samples,
                                        win_params=win_params,
                                        num_wins_to_process=num_wins_to_process)

    pis = {labels.LABELS[pop]: values for pop, values in pis.items()}

    plot_path = config.FIG_PI_DISTRIBUTIONS

    color_schema = colors.ColorSchema(colors.CLASSIFICATION_RANK1_COLORS)

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    plot_distributions_per_pop(pis, axes, color_schema=color_schema)
    axes.legend()
    axes.set_xlabel('pi')
    axes.set_ylabel('Num. windows')
    matplotlib_support.set_axes_background(axes)

    fig.tight_layout()
    fig.savefig(plot_path)
