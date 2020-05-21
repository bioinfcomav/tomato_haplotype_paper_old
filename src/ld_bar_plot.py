
import config

import os

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import numpy

from scipy import interpolate

import statsmodels.api as sm

from variation.variations import VariationsH5

from ld import calc_lds, _get_lds
from passport import get_sample_passports
from pop_building import get_pops


def plot_ld_bars(variations, pops, dist_for_bar_plot, plot_path,
                 dist_for_ld,
                 num_lds_for_loess,
                 max_maf, max_dist,
                 max_num_lds,
                 min_num_samples,
                 cache_dir=None,
                 pops_to_remove=None):

    if pops_to_remove is None:
        pops_to_remove = []

    lds_per_pop = {}
    for pop, samples in pops.items():
        if pop in pops_to_remove:
            continue
        if len(samples) < min_num_samples:
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

    pops = sorted(lds_per_pop.keys())

    bar_edges = numpy.arange(len(pops) + 1)
    bar_width = bar_edges[1] - bar_edges[0]
    bar_x_poss = (bar_edges[1:] + bar_edges[:-1]) / 2
    bar_heights = [lds_per_pop[pop] for pop in pops]

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    axes.bar(bar_x_poss, bar_heights, bar_width)

    axes.set_xticklabels(pops, rotation=45, horizontalalignment='right')
    axes.set_xticks(bar_x_poss)

    axes.set_ylabel(f'LD at {dist_for_ld} bps')

    fig.tight_layout()
    fig.savefig(str(plot_path))

if __name__ == '__main__':

    max_maf = 0.95
    rank = config.RANK1
    max_dist = 1e6
    max_dists_for_decay = [10000, 100000]
    num_lds_for_loess = 5000
    num_lds_for_background = 5000
    max_num_lds = 500000
    min_num_samples_in_pop = 15
    min_lds = 1000
    dist_for_bar_plot = 10000
    cache_dir = config.CACHE_DIR

    sample_passports = get_sample_passports()
    pops_descriptions = {rank: config.ALL_POPS}
    pops = get_pops(pops_descriptions, sample_passports)

    vars_path = config.WORKING_PHASED_H5
    variations = VariationsH5(str(vars_path), 'r')

    pops = {pop:samples for pop, samples in pops.items() if len(samples) > min_num_samples_in_pop}

    out_dir = config.LD_DIR
    os.makedirs(out_dir, exist_ok=True)

    plot_path = out_dir / f'lds_at_{dist_for_bar_plot}_bp.max_maf_{max_maf}.svg'
    plot_ld_bars(variations, pops, dist_for_bar_plot, plot_path, pops_to_remove=[None],
                 dist_for_ld=10000,
                 num_lds_for_loess=num_lds_for_loess,
                 max_maf=max_maf, max_dist=max_dist,
                 min_num_samples=min_num_samples_in_pop,
                 max_num_lds=max_num_lds,
                 cache_dir=cache_dir)
