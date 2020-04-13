
import config

import random
from pprint import pprint
from collections import defaultdict

import numpy

from variation import GT_FIELD, MISSING_INT
from variation.variations import VariationsH5
from variation.variations.filters import VariableAndNotAllMissing, SampleFilter, FLT_VARS

from variation.matrix.stats import counts_and_allels_by_row, counts_by_row

from passport import get_sample_passports
from pop_building import get_pops


class NotEnoughSamplesError(RuntimeError):
    pass


def _set_vars_with_not_enough_gts_to_nan(stat_per_var,
                                         var_has_to_much_missing):
    if issubclass(stat_per_var.dtype.type, numpy.integer):
        stat_per_var = stat_per_var.astype(float)

    stat_per_var[var_has_to_much_missing] = numpy.nan
    return stat_per_var


def calc_pop_stats_per_var(variations, allowed_missing_gts=0,
                           stats_to_calc=None, ploidy=None):
    if stats_to_calc is None:
        stats_to_calc = ['mafs', 'num_alleles', 'unbiased_exp_het',
                         'var_is_poly95', 'var_is_poly75']

    gts = variations[GT_FIELD]

    num_vars = gts.shape[0]
    num_samples = gts.shape[1]
    if ploidy is None:
        ploidy = gts.shape[2]

    allele_counts_per_allele_and_snp, alleles = counts_and_allels_by_row(gts)

    if allele_counts_per_allele_and_snp is None:
        raise ValueError('No genotypes to calc')

    try:
        missing_allele_index = list(alleles).index(MISSING_INT)
    except ValueError:
        missing_allele_index = None

    if missing_allele_index is None:
        num_missing_gts_per_var = numpy.zeros(allele_counts_per_allele_and_snp.shape[0])
    else:
        num_missing_gts_per_var = allele_counts_per_allele_and_snp[:, missing_allele_index]
        allele_counts_per_allele_and_snp = numpy.delete(allele_counts_per_allele_and_snp,
                                                        missing_allele_index, axis=1)

    var_has_to_much_missing = num_missing_gts_per_var > allowed_missing_gts

    if set(stats_to_calc).intersection(['var_is_poly95', 'var_is_poly75', 'mafs', 'unbiased_exp_het']):
        tot_allele_count_per_var = numpy.sum(allele_counts_per_allele_and_snp, axis=1)

    res = {}

    res['num_vars'] = num_vars
    res['num_vars_with_enough_gts'] = num_vars - numpy.sum(var_has_to_much_missing)

    if set(stats_to_calc).intersection(['mafs', 'var_is_poly95', 'var_is_poly75']):
        max_ = numpy.amax(allele_counts_per_allele_and_snp, axis=1)
        # To avoid problems with NaNs
        with numpy.errstate(invalid='ignore'):
            mafs_per_var = max_ / tot_allele_count_per_var
        mafs_per_var = _set_vars_with_not_enough_gts_to_nan(mafs_per_var,
                                                            var_has_to_much_missing)
        if 'mafs' in stats_to_calc:
            res['mafs_per_var'] = mafs_per_var

        if 'var_is_poly95' in stats_to_calc:
            with numpy.errstate(invalid='ignore'):
                res['var_is_poly95'] = mafs_per_var <= 0.95

        if 'var_is_poly75' in stats_to_calc:
            with numpy.errstate(invalid='ignore'):
                res['var_is_poly75'] = mafs_per_var <= 0.75

    if 'num_alleles' in stats_to_calc:
        num_alleles_per_var = numpy.sum(allele_counts_per_allele_and_snp != 0, axis=1)
        res['num_alleles'] = _set_vars_with_not_enough_gts_to_nan(num_alleles_per_var,
                                                                  var_has_to_much_missing)

    if 'unbiased_exp_het' in stats_to_calc:
        with numpy.errstate(invalid='ignore'):
            allele_freq_per_allele_and_var = allele_counts_per_allele_and_snp / tot_allele_count_per_var[:, None]

        exp_het_per_var = 1 - numpy.sum(allele_freq_per_allele_and_var ** ploidy, axis=1)

        num_called_gts_per_var = numpy.full(num_vars, num_samples) - num_missing_gts_per_var
        num_called_gts_per_var_doubled = 2 * num_called_gts_per_var
        unbiased_exp_het = (num_called_gts_per_var_doubled / (num_called_gts_per_var_doubled - 1)) * exp_het_per_var
        res['unbiased_exp_het'] = _set_vars_with_not_enough_gts_to_nan(unbiased_exp_het,
                                                                       var_has_to_much_missing)
    return res


def calc_pop_stats(variations, allowed_missing_gts, percentiles=[25, 50, 75]):
    pop_stas_per_var = calc_pop_stats_per_var(variations, allowed_missing_gts=allowed_missing_gts)
    res = {}
    res['unbiased_exp_het_percentiles'] = numpy.nanpercentile(pop_stas_per_var['unbiased_exp_het'], percentiles)
    res['mean_num_alleles'] = numpy.nanmean(pop_stas_per_var['num_alleles'])
    num_poly95 = numpy.sum(pop_stas_per_var['var_is_poly95'])
    num_poly75 = numpy.sum(pop_stas_per_var['var_is_poly75'])
    res['poly95'] = num_poly95 / pop_stas_per_var['num_vars_with_enough_gts'] * 100
    res['poly75'] = num_poly75 / pop_stas_per_var['num_vars_with_enough_gts'] * 100
    res['ratio_poly75/poly95'] = num_poly75 / num_poly95
    return res


def do_rarefaction_for_population(variations, samples,
                                  rarefaction_range,
                                  allowed_missing_gts=0,
                                  percentiles=[25, 50, 75]):

    if len(samples) < rarefaction_range[0]:
        raise NotEnoughSamplesError()

    max_num_indis = len(samples) + 1 if len(samples) + 1 < rarefaction_range[1] else rarefaction_range[1]
        
    num_indis_range = range(rarefaction_range[0],
                            max_num_indis)

    first = True
    res = {'num_indis': []}
    for num_indis in num_indis_range:
        if num_indis == len(samples):
            samples_for_this_iter = samples
        else:
            samples_for_this_iter = random.sample(samples, num_indis)

        variations_for_this_iter = SampleFilter(samples_for_this_iter)(variations)[FLT_VARS]

        res['num_indis'].append(num_indis)

        pop_stats = calc_pop_stats(variations_for_this_iter,
                                   allowed_missing_gts=allowed_missing_gts,
                                   percentiles=percentiles)
        if first:
            for field in pop_stats.keys():
                res[field] = []
            first = False

        for field, value in pop_stats.items():
            res[field].append(value)

    print(res)
    return res


def keep_variations_variable_in_samples(variations, samples):
    return VariableAndNotAllMissing()(SampleFilter(samples)(variations)[FLT_VARS])[FLT_VARS]


def calc_rarefacted_diversities(variations, pops, rarefaction_range):

    all_samples = {sample for samples in pops.values() for sample in samples}

    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        all_samples = all_samples.intersection(variations.samples)
    all_samples = sorted(all_samples)

    variations = keep_variations_variable_in_samples(variations, all_samples)

    diversities = {}
    for pop, samples in pops.items():

        if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
            samples = sorted(set(samples).intersection(variations.samples))

        try:
            diversities[pop] = do_rarefaction_for_population(variations, samples,
                                                             rarefaction_range=rarefaction_range)
        except NotEnoughSamplesError:
            continue    
    return diversities


if __name__ == '__main__':

    rarefaction_range = (8, 21)

    vars_path = config.WORKING_PHASED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = get_sample_passports()

    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, passports)

    rarefacted_diversities = calc_rarefacted_diversities(variations, pops, rarefaction_range)
    #print(rarefacted_diversities)

    # TODO filtrar a 095 with all_samples
