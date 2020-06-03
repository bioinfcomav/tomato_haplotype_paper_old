
import config

import numpy

from variation import GT_FIELD, MISSING_INT
from variation.variations import VariationsH5
from variation.variations.filters import SampleFilter, FLT_VARS
from variation.matrix.stats import counts_by_row

import passport
import pop_building


def _calc_allele_freq(variations, alleles, min_gt_counts_to_consider_var=10):
    gt_counts = counts_by_row(variations[GT_FIELD], alleles=alleles)

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
                                alleles=None):

    if alleles is None:
        alleles = sorted(set(numpy.unique(variations[GT_FIELD])).difference([MISSING_INT]))

    sample_flt = SampleFilter(samples_in_pop)

    pop_vars = sample_flt(variations)[FLT_VARS]

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


def calc_introgession_freq_for_vars(variations, introgession_config, return_counts=False):
    '''It assumes that target_pop was originated out from a founder_pop and
    that the introgression_source_pop is known
    '''

    samples_in_pop_with_introgressions = introgession_config['samples_in_pop_with_introgressions']
    samples_in_founder_pop = introgession_config['samples_in_founder_pop']
    samples_in_introgression_source_pop = introgession_config['samples_in_introgression_source_pop']

    freq_threshold_to_consider_allele_present_in_founder_pop = introgession_config['freq_threshold_to_consider_allele_present_in_founder_pop']
    freq_threshold_to_consider_allele_common_in_introgression_source_pop = introgession_config['freq_threshold_to_consider_allele_common_in_introgression_source_pop']

    alleles = sorted(set(numpy.unique(variations[GT_FIELD])).difference([MISSING_INT]))

    allele_is_present_in_founder_pop = calc_allele_presence_in_pop(variations, samples_in_founder_pop, freq_threshold_to_consider_allele_present_in_founder_pop, alleles)
    allele_is_not_present_in_founder_pop = numpy.logical_not(allele_is_present_in_founder_pop)
    allele_is_common_in_introgression_origin_pop = calc_allele_presence_in_pop(variations, samples_in_introgression_source_pop, freq_threshold_to_consider_allele_common_in_introgression_source_pop, alleles)

    allele_could_be_introgressed = numpy.logical_and(allele_is_not_present_in_founder_pop,
                                                     allele_is_common_in_introgression_origin_pop)

    samples_in_pop_to_check = samples_in_pop_with_introgressions
    sample_flt = SampleFilter(samples_in_pop_to_check)
    vars_for_pop_to_check = sample_flt(variations)[FLT_VARS]
    res = _calc_allele_freq(vars_for_pop_to_check, alleles)
    allele_freqs = res['freqs']

    if return_counts:
        allele_freqs *= res['max']

    allele_does_not_originate_from_introgression = numpy.logical_not(allele_could_be_introgressed)
    allele_freqs[allele_does_not_originate_from_introgression] = 0
    introgression_freq_per_var = numpy.sum(allele_freqs, axis=1)
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
