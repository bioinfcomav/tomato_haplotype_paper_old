
from numpy.matrixlib.defmatrix import matrix
import config

import numpy
import pandas

from variation import GT_FIELD, MISSING_INT, POS_FIELD
from variation.variations import VariationsH5
from variation.variations.filters import SampleFilter, FLT_VARS
import hashlib, pickle, gzip

import passport
import pop_building
import matrix_methods
import snp_filtering


def _calc_allele_freq(variations, alleles, min_gt_counts_to_consider_var=10):
    #gt_counts = counts_by_row(variations[GT_FIELD], alleles=alleles)
    gt_counts = matrix_methods.counts_and_alleles_by_row(variations[GT_FIELD], alleles=alleles)[0]

    sum_ = numpy.sum(gt_counts, axis=1)
    var_has_not_enough_samples = sum_ < min_gt_counts_to_consider_var
    sum_ = sum_.astype(float)
    sum_[var_has_not_enough_samples] = numpy.nan

    n_alleles = len(alleles)
    sum_ = numpy.repeat(sum_, n_alleles).reshape(gt_counts.shape)

    # To avoid problems with NaNs
    with numpy.errstate(invalid='ignore'):
        freqs = gt_counts / sum_

    return {'freqs': freqs, 'max': numpy.amax(sum_)}


def calc_allele_presence_in_pop(variations, samples_in_pop, freq_threshold_to_consider_allele_present,
                                alleles=None, cache_dir=None):

    if alleles is None:
        alleles = sorted(set(matrix_methods.unique(variations[GT_FIELD])).difference([MISSING_INT]))

    pop_vars = snp_filtering.filter_variations(variations,
                                               samples_to_keep=samples_in_pop,
                                               cache_dir=cache_dir)

    allele_freqs = _calc_allele_freq(pop_vars, alleles)['freqs']

    allele_freqs[numpy.isnan(allele_freqs)] = 0
    allele_is_in_pop = allele_freqs >= freq_threshold_to_consider_allele_present

    return allele_is_in_pop


def _calc_introgession_freq_for_vars(variations, samples_in_pop, allele_could_be_introgressed, alleles, return_counts=False):

    sample_flt = SampleFilter(samples_in_pop)

    pop_vars = sample_flt(variations)[FLT_VARS]

    res = _calc_allele_freq(pop_vars, alleles)
    allele_freqs = res['freqs']

    if return_counts:
        allele_freqs *= res['max']

    allele_could_not_be_introgressed = numpy.logical_not(allele_could_be_introgressed)
    allele_freqs[allele_could_not_be_introgressed] = 0
    introgression_freq_per_var = numpy.sum(allele_freqs, axis=1)
    return introgression_freq_per_var


def _dict_to_str(dict_):
    str_ = ''
    for key in sorted(dict_, key=str):
        str_ += str(dict_[key])
    return str_


def calc_introgession_freq_for_vars(variations, introgession_config,
                                    return_counts=False,
                                    cache_dir=None):
    '''It assumes that target_pop was originated out from a founder_pop and
    that the introgression_source_pop is known
    '''

    cache_path = None
    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += str(variations.num_variations)
        key += 'introgession_config' + _dict_to_str(introgession_config)
        key += 'return_counts' + str(return_counts)
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('introgression_freqs_for_vars.' + key + '.pickle')
        if cache_path.exists():
            return pickle.load(gzip.open(cache_path, 'rb'))

    samples_in_pop_with_introgressions = introgession_config['samples_in_pop_with_introgressions']
    samples_in_founder_pop = introgession_config['samples_in_founder_pop']
    samples_in_introgression_source_pop = introgession_config['samples_in_introgression_source_pop']

    freq_threshold_to_consider_allele_present_in_founder_pop = introgession_config['freq_threshold_to_consider_allele_present_in_founder_pop']
    freq_threshold_to_consider_allele_common_in_introgression_source_pop = introgession_config['freq_threshold_to_consider_allele_common_in_introgression_source_pop']

    alleles = sorted(set(matrix_methods.unique(variations[GT_FIELD])).difference([MISSING_INT]))

    allele_is_present_in_founder_pop = calc_allele_presence_in_pop(variations,
                                                                   samples_in_founder_pop,
                                                                   freq_threshold_to_consider_allele_present_in_founder_pop,
                                                                   alleles,
                                                                   cache_dir=cache_dir)
    allele_is_not_present_in_founder_pop = numpy.logical_not(allele_is_present_in_founder_pop)
    allele_is_common_in_introgression_origin_pop = calc_allele_presence_in_pop(variations,
                                                                               samples_in_introgression_source_pop,
                                                                               freq_threshold_to_consider_allele_common_in_introgression_source_pop,
                                                                               alleles,
                                                                               cache_dir=cache_dir)

    allele_could_be_introgressed = numpy.logical_and(allele_is_not_present_in_founder_pop,
                                                     allele_is_common_in_introgression_origin_pop)

    allele_does_not_originate_from_introgression = numpy.logical_not(allele_could_be_introgressed)

    samples_in_pop_to_check = samples_in_pop_with_introgressions
    vars_for_pop_to_check = snp_filtering.filter_variations(variations,
                                                            samples_to_keep=samples_in_pop_to_check,
                                                            cache_dir=cache_dir)

    res = _calc_allele_freq(vars_for_pop_to_check, alleles)
    allele_freqs = res['freqs']

    if return_counts:
        allele_freqs *= res['max']

    allele_freqs[allele_does_not_originate_from_introgression] = 0
    introgression_freq_per_var = numpy.sum(allele_freqs, axis=1)

    assert variations.num_variations == introgression_freq_per_var.size
    introgression_freq_per_var = pandas.Series(introgression_freq_per_var,
                                               index=numpy.arange(introgression_freq_per_var.size))

    if cache_dir:
        pickle.dump(introgression_freq_per_var, gzip.open(cache_path, 'wb'))

    return introgression_freq_per_var


if __name__ == '__main__':
    vars_path = config.WORKING_H5
    variations = VariationsH5(str(vars_path), 'r')

    target_pop = 'slc_pe'
    founder_pop = 'slc_ma'
    introgression_source_pop = 'sp_pe'

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, passports)

    introgession_config = {'samples_in_pop_with_introgressions': pops[target_pop],
                           'samples_in_founder_pop': pops[founder_pop],
                           'samples_in_introgression_source_pop': pops[introgression_source_pop],
                           'freq_threshold_to_consider_allele_present_in_founder_pop' : config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_PRESENT,
                           'freq_threshold_to_consider_allele_common_in_introgression_source_pop': config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_COMMON,
                          }
    calc_introgession_freq_for_vars(variations, introgession_config)
