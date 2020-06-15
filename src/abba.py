
import config

from collections import OrderedDict
import itertools
import math
import random
from array import array

import numpy

from scipy.stats import norm, zscore

from variation import GT_FIELD, MISSING_INT
from variation.variations import VariationsH5, VariationsArrays
from variation.variations.filters import NonBiallelicFilter, FLT_VARS, SampleFilter
from variation.variations.pipeline import Pipeline
from variation.matrix.stats import counts_by_row

import passport
import pop_building


def _filter_vars(variations, samples_to_keep):

    chunk_size = 10000
    max_chunks_to_process = None

    biallelic_filter = NonBiallelicFilter(samples=samples_to_keep)

    pipeline = Pipeline()
    pipeline.append(biallelic_filter)
    filtered_vars = VariationsArrays()

    result = pipeline.run(variations, filtered_vars, chunk_size=chunk_size,
                            max_chunks_to_process=max_chunks_to_process,
                            kept_fields=[GT_FIELD])
    return filtered_vars


def calc_ab_freq_products_for_patterns(variations, pops_in_ladderized_order, patterns=None):

    if patterns is None:
        patterns = [tuple(pattern) for pattern in itertools.product('AB', repeat=len(pops_in_ladderized_order)) if len(set(pattern)) > 1]

    sorted_pop_names = list(pops_in_ladderized_order.keys())
    origin_pop = sorted_pop_names[-1]

    samples = [sample for samples in pops_in_ladderized_order.values() for sample in samples]

    variations = _filter_vars(variations, samples)

    alleles = sorted([alelle for alelle in numpy.unique(variations[GT_FIELD].flat) if alelle != MISSING_INT])

    vars_for_pops = {pop: SampleFilter(samples)(variations)[FLT_VARS] for pop, samples in pops_in_ladderized_order.items()}

    allele_freqs_for_pops = {}
    for pop, vars_for_this_pop in vars_for_pops.items():
        counts = counts_by_row(vars_for_this_pop[GT_FIELD], alleles=alleles)
        freqs = counts / counts.sum(axis=1)[:, numpy.newaxis]
        allele_freqs_for_pops[pop] = freqs

    a_allele_index = numpy.argmax(allele_freqs_for_pops[origin_pop], axis=1)

    snp_indexes = numpy.arange(a_allele_index.size)

    freqs = {'A': {}, 'B': {}}
    for pop, allele_freqs in allele_freqs_for_pops.items():
        a_freq = allele_freqs[snp_indexes, a_allele_index]
        freqs['A'][pop] = a_freq
        freqs['B'][pop] = 1 - a_freq

    products = {}
    for pattern in patterns:
        product = None
        for allele, pop in zip(pattern, sorted_pop_names):
            if product is None:
                product = freqs[allele][pop]
            else:
                product *= freqs[allele][pop]
        products[pattern] = product
    return products


def calc_d_statistic(variations, pops_in_ladderized_order, num_jacknife_repeats=None,
                     jackknife_rate_to_remove=0.3):
    if len(pops_in_ladderized_order) != 4:
        raise ValueError('D statistic is defined for 4 pops')

    patterns = [('A', 'B', 'B', 'A'), ('B', 'A', 'B', 'A')]
    products = calc_ab_freq_products_for_patterns(variations, pops_in_ladderized_order, patterns=patterns)
    abba = products[('A', 'B', 'B', 'A')]
    baba = products[('B', 'A', 'B', 'A')]

    bottom_per_snp = abba + baba
    top_per_snp = abba - baba

    bottom = numpy.sum(bottom_per_snp)
    top = numpy.sum(top_per_snp)

    if math.isclose(bottom, 0):
        d_stat = 0
    else:
        d_stat = top / bottom
    result = {'d_stat': d_stat}

    if num_jacknife_repeats is not None:
        num_snps = abba.size
        last_posible_start_point = int(num_snps * (1- jackknife_rate_to_remove))
        d_stats = array('f')
        for _ in range(num_jacknife_repeats):
            start = random.randint(0, last_posible_start_point)
            end = start + last_posible_start_point
            mask1 = slice(start)
            mask2 = slice(end, None)
            this_bottom_per_snp = numpy.concatenate([bottom_per_snp[mask1], bottom_per_snp[mask2]])
            this_top_per_snp = numpy.concatenate([top_per_snp[mask1], top_per_snp[mask2]])
            this_bottom = numpy.sum(this_bottom_per_snp)
            if math.isclose(bottom, 0):
                this_d_stat = 0
            else:
                this_d_stat = numpy.sum(this_top_per_snp) / this_bottom
            d_stats.append(this_d_stat)
        d_stats = numpy.array(d_stats)
        result['jackknifed_d_stats'] = d_stats

        # p value calculated according to
        # http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics/
        d_standard_error = numpy.std(d_stats) / math.sqrt(d_stats.size)
        d_z_score = numpy.mean(d_stats) / d_standard_error
        d_p_value = 2 * norm.sf(abs(d_z_score))
        result['d_p_value'] = d_p_value

    return result


if __name__ == '__main__':

    num_jacknife_repeats = 1000

    vars_path = config.WORKING_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()
    pops = pop_building.get_pops({config.RANK1: config.ALL_POPS}, passports)

    pops_names_for_tree = [['slc_pe', 'slc_ma', 'sp_ec', 'sp_pe'],
                           ['slc_ec', 'slc_ma', 'sp_ec', 'sp_pe'],
                           ['slc_ma', 'slc_ec', 'sp_ec', 'sp_pe'],
                           ['slc_ma', 'sp_ec', 'slc_ec', 'sp_pe'],
                           ['slc_ec', 'slc_pe', 'slc_ma', 'sp_ec'],
                           ['sp_ec', 'slc_ec', 'slc_ma', 'sp_pe']]

    for pop_names in pops_names_for_tree:
        pops_for_tree = OrderedDict([(pop, pops[pop]) for pop in pop_names])
        res = calc_d_statistic(variations, pops_for_tree, num_jacknife_repeats=num_jacknife_repeats)
        print(pop_names, res['d_stat'], res['d_p_value'])

    pops_for_tree = ['slc_ec', 'sll_mx', 'slc_pe', 'slc_ma', 'sp_ec', 'sp_pe_inter-andean', 'sp_pe']
    pops_for_tree = OrderedDict([(pop, pops[pop]) for pop in pops_for_tree])
    products = calc_ab_freq_products_for_patterns(variations, pops_for_tree)

    sum_ = None
    for product in products.values():
        if sum_ is None:
            sum_ = product
        else:
            sum_ += product
    print(numpy.nansum(sum_))
