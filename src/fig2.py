
import config

from pprint import pprint
import hashlib, pickle, gzip

import numpy
from scipy import interpolate

import statsmodels.api as sm

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib

from variation.variations import VariationsH5

from pop_building import get_pops
from passport import get_sample_passports
from diversities_vars import calc_var_diversities as calc_var_diversities_orig
from diversities_haplos import calc_haplo_diversities as calc_haplo_diversities_orig
from ld import _get_lds
import colors
import labels


Y_LABEL_SIZE = 15
X_TICK_LABEL_SIZE = 12


def calc_var_diversities(variations, pops, num_samples, num_repeats, cache_dir=None):
    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += 'num_variations' + str(variations.num_variations)
        key += 'pops' + str(sorted(pops.items()))
        key += 'num_samples' + str(num_samples)
        key += 'num_repeats' + str(num_repeats)
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('var_diversities' + key + '.pickle.gz')
        if cache_path.exists():
            return pickle.load(gzip.open(cache_path, 'rb'))

    res = calc_var_diversities_orig(variations=variations, pops=pops, num_samples=num_samples, num_repeats=num_repeats)

    if cache_dir:
        pickle.dump(res, gzip.open(cache_path, 'wb'))
    return res


def calc_haplo_diversities(variations, pops, num_samples, num_repeats,
                           win_params, cache_dir=None):
    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += 'num_variations' + str(variations.num_variations)
        key += 'pops' + str(sorted(pops.items()))
        key += 'num_samples' + str(num_samples)
        key += 'num_repeats' + str(num_repeats)
        key += 'win_params' + str(win_params)
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('haplo_diversities' + key + '.pickle.gz')
        if cache_path.exists():
            return pickle.load(gzip.open(cache_path, 'rb'))
    
    res = calc_haplo_diversities_orig(imputed_variations, pops=pops, num_samples=num_samples, num_repeats=num_repeats,
                                      win_params=win_params)
    if cache_dir:
        pickle.dump(res, gzip.open(cache_path, 'wb'))
    return res


def calc_ld(variations, pops, dist_for_ld,
            num_lds_for_loess,
            max_maf, max_dist,
            max_num_lds,
            min_num_samples_in_pop,
            cache_dir=None,
            pops_to_remove=None):

    if pops_to_remove is None:
        pops_to_remove = []

    lds_per_pop = {}
    for pop, samples in pops.items():
        if pop in pops_to_remove:
            continue
        if len(samples) < min_num_samples_in_pop:
            continue

        res = _get_lds(variations, samples, max_maf=max_maf,
                        max_dist=max_dist,
                        min_num_samples=min_num_samples_in_pop,
                        max_num_lds=max_num_lds,
                        cache_dir=cache_dir)
        lds = res['lds']
        dists = res['dists']
        lowess_dists = dists[:num_lds_for_loess]
        lowess_lds = lds[:num_lds_for_loess]
        delta = 0.01 * (max(lowess_lds) - min(lowess_dists))
        lowess = sm.nonparametric.lowess(lowess_lds, lowess_dists, delta=delta)
        lowess_dists, lowess_lds = list(zip(*lowess))
        ld_funct = interpolate.interp1d(lowess_dists, lowess_lds)
        try:
            ld = ld_funct([dist_for_ld])[0]
        except ValueError:
            continue
        lds_per_pop[pop] = ld
    return lds_per_pop


def _plot_diversities(diversity_dict, pop_order, axes, color_schema=None):

    heights = [diversity_dict[pop]['mean'] for pop in pop_order]
    x_values = numpy.arange(len(pop_order))

    yerr = [diversity_dict[pop]['cim'] for pop in pop_order]

    if color_schema is None:
        color = None
    else:
        color = [color_schema[pop] for pop in pop_order]

    axes.bar(x_values, heights, yerr=yerr, color=color)
    return {'x_values': x_values, 'x_labels': pop_order}


def turn_off_x_axis(axes):
    axes.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=False)


def set_y_ticks_right(axes):
    axes.tick_params(axis='y', which='both', right=True,
                     left=False, labelleft=False, labelright=True)
    axes.yaxis.set_label_position('right')


def set_axis_color(axes, color='black'):
    for child in axes.get_children():
        if isinstance(child, matplotlib.spines.Spine):
            child.set_color('#dddddd')


def  plot_fig2(var_diversities, haplo_diversities, lds, plot_path):

    pop_order = ['sp_pe', 'sp_x_sp', 'sp_pe_inter-andean',
                           'sp_ec',
                           'slc_co',
                           'slc_ma',
                           'slc_world',
                           'sll_mx', 'sll_vint', 'sll_old_cultivars',
                           'sll_modern',
                           'sp_x_sl',
                           None,
                           'slc_pe',
                           'slc_ec']

    color_schema = colors.ColorSchema(colors.CLASSIFICATION_RANK1_COLORS)

    fig = Figure((10, 10))
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes1 = fig.add_subplot(221)

    poly95 = var_diversities['num_poly95']
    pops = sorted(poly95.keys(), key=lambda x: pop_order.index(x))

    _plot_diversities(poly95, pops, axes1, color_schema=color_schema)
    turn_off_x_axis(axes1)
    axes1.set_ylabel('Num. poly variations (95%)', fontsize=Y_LABEL_SIZE)

    axes2 = fig.add_subplot(222)
    _plot_diversities(haplo_diversities['num_uniq_haplos'], pops, axes2, color_schema=color_schema)
    turn_off_x_axis(axes2)
    set_y_ticks_right(axes2)
    axes2.set_ylabel('Mean num. unique haplotypes', fontsize=Y_LABEL_SIZE)

    axes3 = fig.add_subplot(224)
    res = _plot_diversities(haplo_diversities['num_variable_vars'], pops, axes3, color_schema=color_schema)
    set_y_ticks_right(axes3)
    axes3.set_ylabel('Mean num. variable variations in genome region', fontsize=Y_LABEL_SIZE)

    axes3.set_xticklabels(res['x_labels'], rotation=45, ha='right', fontsize=X_TICK_LABEL_SIZE)
    axes3.set_xticks([labels.get_long_label(pop) for pop in res['x_values']])

    axes4 = fig.add_subplot(223)
    heights = [lds.get(pop, 0) for pop in pops]
    x_values = numpy.arange(len(pops))
    color = [color_schema[pop] for pop in pops]
    axes4.bar(x_values, heights, color=color)
    axes4.set_ylabel('LD at 10 Kb', fontsize=Y_LABEL_SIZE)
    axes4.set_xticklabels([labels.get_long_label(pop) for pop in pops], rotation=45, ha='right', fontsize=X_TICK_LABEL_SIZE)
    axes4.set_xticks(x_values)

    axes1.set_facecolor('white')
    axes2.set_facecolor('white')
    axes3.set_facecolor('white')
    axes4.set_facecolor('white')
    set_axis_color(axes1)
    set_axis_color(axes2)
    set_axis_color(axes3)
    set_axis_color(axes4)

    fig.tight_layout()
    fig.savefig(str(plot_path))


if __name__  == '__main__':

    out_dir = config.PAPER_FIGURES_DIR
    out_dir.mkdir(exist_ok=True)
    plot_path = out_dir / 'fig2.svg'

    num_repeats = 100
    percent_indiv_used = 0.75
    min_num_indi_in_pop = 10

    pops_to_exclude = [None]

    dist_for_ld = 10000
    num_lds_for_loess = 5000
    max_maf_for_ld = 0.95
    max_dist_for_ld = 1e6
    max_num_lds = 500000

    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}

    cache_dir = config.CACHE_DIR

    vars_path = config.WORKING_H5
    variations = VariationsH5(str(vars_path), 'r')

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    imputed_variations = VariationsH5(str(vars_path), 'r')

    passports = get_sample_passports()

    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, passports)

    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        pops = {pop:sorted(set(samples).intersection(variations.samples))  for pop, samples in pops.items()}

    pops = {pop: sorted(samples) for pop, samples in pops.items() if len(samples) > min_num_indi_in_pop and pop not in pops_to_exclude}

    num_samples_per_pop ={pop: len(samples) for pop, samples in pops.items()}

    num_samples = round(min(num_samples_per_pop.values()) * percent_indiv_used)

    var_diversities = calc_var_diversities(variations, pops, num_samples, num_repeats=num_repeats, cache_dir=cache_dir)

    if True:
        haplo_diversities = calc_haplo_diversities(imputed_variations, pops, num_samples, num_repeats=num_repeats,
                                                   win_params=win_params, cache_dir=cache_dir)
    else:
        haplo_diversities = None

    lds = calc_ld(variations, pops=pops, dist_for_ld=dist_for_ld,
                  max_maf=max_maf_for_ld,
                  num_lds_for_loess=num_lds_for_loess,
                  max_dist=max_dist_for_ld,
                  max_num_lds=max_num_lds,
                  min_num_samples_in_pop=num_samples,
                  cache_dir=cache_dir)

    plot_fig2(var_diversities, haplo_diversities, lds, plot_path)

