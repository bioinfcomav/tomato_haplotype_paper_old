
import config

import itertools
import os
import gzip, pickle, hashlib
import math

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import statsmodels.api as sm

from variation.variations import VariationsH5

from passport import get_sample_passports
from pop_building import get_pops
from ld import calc_lds, _get_lds


def plot_ld_decay(lds, dists, plot_path,
                  max_dist_for_decay=None, num_lds_for_loess=None,
                  num_lds_for_background=None,
                  min_lds=None):

    if max_dist_for_decay:
        mask = dists < max_dist_for_decay
        dists = dists[mask]
        lds = lds[mask]

    if min_lds is not None:
        if lds.size < min_lds:
            raise RuntimeError('Not enough lds')

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)
    axes.scatter(dists[:num_lds_for_background], lds[:num_lds_for_background], alpha=0.3)

    lowess_dists = dists[:num_lds_for_loess]
    lowess_lds = lds[:num_lds_for_loess]
    delta = 0.01 * (max(lowess_lds) - min(lowess_dists))
    lowess = sm.nonparametric.lowess(lowess_lds, lowess_dists, delta=delta)
    lowess_x, lowess_y = list(zip(*lowess))

    axes.plot(lowess_x, lowess_y)

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

    for max_dist_for_decay in max_dists_for_decay:
        pop_dir = out_dir / 'ld_for_pops'
        os.makedirs(pop_dir, exist_ok=True)

        for pop, samples in pops.items():
            res = _get_lds(variations, samples, max_maf=max_maf,
                           max_dist=max_dist,
                           min_num_samples=min_num_samples_in_pop,
                           max_num_lds=max_num_lds,
                           cache_dir=cache_dir)
            plot_path = pop_dir / f'ld_{pop}.svg'
            print(plot_path)
            try:
                plot_ld_decay(res['lds'], res['dists'], plot_path,
                              max_dist_for_decay=max_dist_for_decay,
                              num_lds_for_loess=num_lds_for_loess,
                              num_lds_for_background=num_lds_for_background,
                              min_lds=min_lds)
            except RuntimeError:
                continue
