
import config

import random
from pprint import pprint
from collections import defaultdict
import os

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation import GT_FIELD, MISSING_INT
from variation.variations import VariationsH5
from variation.variations.filters import VariableAndNotAllMissing, SampleFilter, FLT_VARS

from variation.matrix.stats import counts_and_allels_by_row, counts_by_row

from passport import get_sample_passports
from pop_building import get_pops
import colors
from snp_filtering import keep_variations_variable_in_samples


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
                         'var_is_poly95', 'var_is_poly80',
                         'var_is_variable']

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

    if set(stats_to_calc).intersection(['var_is_poly95', 'var_is_poly80', 'mafs', 'unbiased_exp_het']):
        tot_allele_count_per_var = numpy.sum(allele_counts_per_allele_and_snp, axis=1)

    res = {}

    res['num_vars'] = num_vars
    res['num_vars_with_enough_gts'] = num_vars - numpy.sum(var_has_to_much_missing)

    if set(stats_to_calc).intersection(['mafs', 'var_is_poly95', 'var_is_poly80', 'var_is_variable']):
        max_ = numpy.amax(allele_counts_per_allele_and_snp, axis=1)
        # To avoid problems with NaNs
        with numpy.errstate(invalid='ignore'):
            mafs_per_var = max_ / tot_allele_count_per_var
        mafs_per_var = _set_vars_with_not_enough_gts_to_nan(mafs_per_var,
                                                            var_has_to_much_missing)
        if 'mafs' in stats_to_calc:
            res['mafs_per_var'] = mafs_per_var

        if 'var_is_variable' in stats_to_calc:
            var_has_enough_data = numpy.logical_not(var_has_to_much_missing)
            var_is_variable = numpy.logical_not(numpy.isclose(mafs_per_var, 1))
            var_is_variable = numpy.logical_and(var_is_variable, var_has_enough_data)
            res['var_is_variable'] = var_is_variable

        if 'var_is_poly95' in stats_to_calc:
            with numpy.errstate(invalid='ignore'):
                res['var_is_poly95'] = mafs_per_var <= 0.95

        if 'var_is_poly80' in stats_to_calc:
            with numpy.errstate(invalid='ignore'):
                res['var_is_poly80'] = mafs_per_var <= 0.75

    if 'num_alleles' in stats_to_calc:
        num_alleles_per_var = numpy.sum(allele_counts_per_allele_and_snp != 0, axis=1)
        res['num_alleles'] = _set_vars_with_not_enough_gts_to_nan(num_alleles_per_var,
                                                                  var_has_to_much_missing)

    if 'unbiased_exp_het' in stats_to_calc:
        num_called_gts_per_var = numpy.full(num_vars, num_samples) - num_missing_gts_per_var
        with numpy.errstate(invalid='ignore'):
            allele_freq_per_allele_and_var = allele_counts_per_allele_and_snp / tot_allele_count_per_var[:, None]

        exp_het_per_var = 1 - numpy.sum(allele_freq_per_allele_and_var ** ploidy, axis=1)

        num_called_gts_per_var_doubled = 2 * num_called_gts_per_var
        unbiased_exp_het = (num_called_gts_per_var_doubled / (num_called_gts_per_var_doubled - 1)) * exp_het_per_var
        res['unbiased_exp_het'] = _set_vars_with_not_enough_gts_to_nan(unbiased_exp_het,
                                                                       var_has_to_much_missing)
    res['stats_calc'] = stats_to_calc

    return res


def calc_pop_stats(variations, allowed_missing_gts=0, percentiles=[25, 50, 75], ploidy=None,
                   stats_to_calc=None):
    if stats_to_calc is None:
        stats_to_calc_per_var = None
    elif not set(stats_to_calc).difference(['unbiased_exp_het_mean']):
        stats_to_calc_per_var = ['unbiased_exp_het']
    else:
        raise NotImplementedError()

    pop_stas_per_var = calc_pop_stats_per_var(variations,
                                              allowed_missing_gts=allowed_missing_gts,
                                              ploidy=ploidy)
    res = {}
    if 'unbiased_exp_het' in pop_stas_per_var['stats_calc']:
        res['unbiased_exp_het_percentiles'] = numpy.nanpercentile(pop_stas_per_var['unbiased_exp_het'], percentiles)
    if 'unbiased_exp_het' in pop_stas_per_var['stats_calc']:
        res['unbiased_exp_het_mean'] = numpy.nanmean(pop_stas_per_var['unbiased_exp_het'])
    if 'num_alleles' in pop_stas_per_var['stats_calc']:
        res['mean_num_alleles'] = numpy.nanmean(pop_stas_per_var['num_alleles'])
    if 'var_is_variable' in pop_stas_per_var['stats_calc']:
        res['num_variable_vars'] = numpy.sum(pop_stas_per_var['var_is_variable'])
    if 'var_is_poly95' in pop_stas_per_var['stats_calc']:
        res['num_poly95'] = numpy.sum(pop_stas_per_var['var_is_poly95'])
        if 'num_vars_with_enough_gts' in pop_stas_per_var:
            res['poly95'] = res['num_poly95'] / pop_stas_per_var['num_vars_with_enough_gts'] * 100
    if 'var_is_poly80' in pop_stas_per_var['stats_calc']:
        res['num_poly80'] = numpy.sum(pop_stas_per_var['var_is_poly80'])
        if 'num_vars_with_enough_gts' in pop_stas_per_var:
            res['poly80'] = res['num_poly80'] / pop_stas_per_var['num_vars_with_enough_gts'] * 100
    try:
        res['ratio_poly80/poly95'] = res['num_poly80'] / res['num_poly95']
    except KeyError:
        pass
    try:
        res ['pi'] = numpy.sum(pop_stas_per_var['unbiased_exp_het']) / pop_stas_per_var['num_vars_with_enough_gts']
    except KeyError:
        pass
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
    res = {'num_samples': []}
    for num_indis in num_indis_range:
        if num_indis == len(samples):
            samples_for_this_iter = samples
        else:
            samples_for_this_iter = random.sample(samples, num_indis)

        variations_for_this_iter = SampleFilter(samples_for_this_iter)(variations)[FLT_VARS]

        res['num_samples'].append(num_indis)

        pop_stats = calc_pop_stats(variations_for_this_iter,
                                   allowed_missing_gts=allowed_missing_gts,
                                   percentiles=percentiles)
        if first:
            for field in pop_stats.keys():
                res[field] = []
            first = False

        for field, value in pop_stats.items():
            res[field].append(value)

    return res


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


def _plot_rarefacted_diversities_for_pops(diversities_per_pop, num_samples_per_pop, plot_path, pop_colors=None, y_lims=None):
    pops = sorted(diversities_per_pop.keys(), key=str)

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    color_schema = colors.ColorSchema(pop_colors)

    for pop in pops:
        color = color_schema[pop]
        x_values = num_samples_per_pop[pop]
        pop_values_for_diversity_field = diversities_per_pop[pop]

        if isinstance(pop_values_for_diversity_field[0], (float, int)):
            y_values = pop_values_for_diversity_field
            axes.plot(x_values, y_values, label=pop, color=color)
        elif isinstance(pop_values_for_diversity_field[0], numpy.ndarray) and pop_values_for_diversity_field[0].size == 3:
            quartiles_low = [quartile_values[0] for quartile_values in pop_values_for_diversity_field]
            medians = [quartile_values[1] for quartile_values in pop_values_for_diversity_field]
            quartiles_high = [quartile_values[-1] for quartile_values in pop_values_for_diversity_field]
            axes.plot(x_values, medians, label=pop, color=color, zorder=10)
            axes.fill_between(x_values, quartiles_low, quartiles_high, alpha=0.5, color=color, zorder=5)

    if y_lims is not None:
        axes.set_ylim(*y_lims)

    axes.legend()

    fig.tight_layout()
    fig.savefig(str(plot_path))


def plot_rarefacted_diversities(rarefacted_diversities, out_dir, pop_colors=None, only_this_pops=None, y_lims=None):
    if only_this_pops is None:
        pops = sorted(rarefacted_diversities.keys(), key=str)
    else:
        pops = only_this_pops

    if y_lims is None:
        y_lims = {}

    diversity_fields = {field for pop_diversities in rarefacted_diversities.values() for field in pop_diversities.keys()}
    diversity_fields = diversity_fields.difference(['num_samples'])

    for field in diversity_fields:
        diversities_per_pop = {pop: rarefacted_diversities[pop][field] for pop in pops}
        num_samples_per_pop = {pop: rarefacted_diversities[pop]['num_samples'] for pop in pops}
        field_no_slash = field.replace('/', '_')

        pops_str = '-'.join(map(str, pops))
        plot_path = out_dir / f'{field_no_slash}.{pops_str}.svg'

        _plot_rarefacted_diversities_for_pops(diversities_per_pop, num_samples_per_pop, plot_path, pop_colors=pop_colors, y_lims=y_lims.get(field))


if __name__ == '__main__':

    rarefaction_range = (8, 23)

    vars_path = config.WORKING_PHASED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = get_sample_passports()

    main_pops = ['sp_pe' ,'sp_ec', 'slc_ma', 'slc_ec', 'slc_pe', 'sll_mx']
    vintage_pops = ['sll_old_cultivars', 'sll_vint', 'slc_world']
    hybrid_pops = ['sll_modern', 'sp_x_sl', 'sp_x_sp']
    all_pops = main_pops + vintage_pops + hybrid_pops

    pops_descriptions = {config.RANK1: all_pops}
    pops = get_pops(pops_descriptions, passports)
    pop_colors = colors.CLASSIFICATION_RANK1_COLORS

    rarefacted_diversities = calc_rarefacted_diversities(variations, pops, rarefaction_range)
    
    out_dir = config.RAREFACTION_VARS_DIR
    os.makedirs(out_dir, exist_ok=True)

    y_lims = {'mean_num_alleles': (1, 2),
              'poly80': (0, 20),
              'poly95': (0, 70),
              'ratio_poly80/poly95': (0, 0.6),
              'unbiased_exp_het_percentiles': (0, 0.35)
              }

    plot_rarefacted_diversities(rarefacted_diversities, out_dir, pop_colors, only_this_pops=main_pops, y_lims=y_lims)
    plot_rarefacted_diversities(rarefacted_diversities, out_dir, pop_colors, only_this_pops=vintage_pops, y_lims=y_lims)
    plot_rarefacted_diversities(rarefacted_diversities, out_dir, pop_colors, only_this_pops=hybrid_pops, y_lims=y_lims)

    # TODO filtrar a 095 with all_samples?
