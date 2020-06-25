
import config

from collections import Counter

import numpy
import pandas

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import single, dendrogram

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation import GT_FIELD, MISSING_INT
from variation.variations import VariationsH5
from variation.variations.filters import SampleFilter, FLT_VARS

import passport
import pop_building
import matplotlib_support


def _allele_freqs_to_dframe(freqs):
    pops = sorted(freqs.keys(), key=str)
    alleles = sorted({allele for pop_freqs in freqs.values() for allele in pop_freqs.keys()})

    freq_array = numpy.zeros((len(alleles), len(pops)))
    freq_dframe = pandas.DataFrame(freq_array, index=alleles, columns=pops)

    for pop, pop_freqs in freqs.items():
        for allele, freq in pop_freqs.items():
            freq_dframe.loc[allele, pop] = freq

    assert numpy.allclose(freq_dframe.sum(axis=0).values, 1)
    return freq_dframe


def calc_allele_freq(variations, pops):

    if MISSING_INT in variations[GT_FIELD]:
        raise ValueError('Allele freq calculation not implemented with missing gts')

    freqs = {}
    for pop, samples in pops.items():
        vars_for_this_pop = SampleFilter(samples)(variations)[FLT_VARS]
        gts = vars_for_this_pop[GT_FIELD]
        counts = Counter()
        for ploid_id in range(gts.shape[2]):
            alleles, this_counts = numpy.unique(gts[:, :, ploid_id], axis=1, return_counts=True)
            for allele_idx, count in enumerate(this_counts):
                allele = tuple(alleles[:, allele_idx])
                counts[allele] += count
    
        pop_freqs = {allele: count / sum(counts.values()) for allele, count in counts.items()}
        freqs[pop] = pop_freqs

    freqs = _allele_freqs_to_dframe(freqs)

    return freqs


def sort_dframe_columns_by_similarity(dframe):

    dist = pdist(dframe.T)
    tree = single(dist)
    leaves_order = dendrogram(tree, no_plot=True)['leaves']

    columns = list(dframe.columns)
    new_columns = [columns[idx] for idx in leaves_order]

    return dframe.reindex(columns=new_columns)


def plot_allele_freqs_bars(freqs, axes, pop_order=None, col_width=0.9):

    if pop_order is None:
        freqs = sort_dframe_columns_by_similarity(freqs)
    sorted_pops = list(freqs.columns)

    alleles = list(freqs.sum(axis=1).sort_values(ascending=False).index)
    orig_allele_order = list(freqs.index)

    num_pops = freqs.shape[1]
    x_values = numpy.arange(num_pops)
    width = (x_values[1] - x_values[0]) * col_width

    bottom = None
    for allele in alleles:
        allele_idx = orig_allele_order.index(allele)

        heights = freqs.iloc[allele_idx, :].values
        axes.bar(x_values, heights, bottom=bottom, width=width, label=''.join(map(str, allele)))

        if bottom is None:
            bottom = heights
        else:
            bottom += heights

    matplotlib_support.set_x_ticks(x_values, sorted_pops, axes=axes, rotation=45)

    return {'alleles': alleles, 'x_values': x_values, 'sorted_pops': sorted_pops}


if __name__ == '__main__':

    k_range = (2, 11)
    max_maf = 0.95

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    variations = variations.get_chunk(slice(2000, 2005))

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, passports)

    freqs = calc_allele_freq(variations, pops)

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    plot_allele_freqs_bars(freqs, axes)

    plot_path = config.TMP_DIR / 'allele_freqs.svg'
    fig.tight_layout()
    fig.savefig(plot_path)
