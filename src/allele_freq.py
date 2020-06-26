
import config

from collections import Counter
import random

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


def _count_alleles(gts):
    counts = Counter()
    for ploid_id in range(gts.shape[2]):
        alleles, this_counts = numpy.unique(gts[:, :, ploid_id], axis=1, return_counts=True)
        for allele_idx, count in enumerate(this_counts):
            allele = tuple(alleles[:, allele_idx])
            counts[allele] += count
    return counts


def calc_allele_freq(variations, pops):

    if MISSING_INT in variations[GT_FIELD][:]:
        raise ValueError('Allele freq calculation not implemented with missing gts')

    freqs = {}
    for pop, samples in pops.items():
        vars_for_this_pop = SampleFilter(samples)(variations)[FLT_VARS]
        gts = vars_for_this_pop[GT_FIELD]
        counts = _count_alleles(gts)
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


def plot_allele_freqs_bars(freqs, axes, pop_order=None, col_width=0.9, color_schema=None):

    if pop_order is None:
        freqs = sort_dframe_columns_by_similarity(freqs)
    else:
        freqs = freqs.reindex(columns=pop_order)
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

        if color_schema is None:
            color = None
        else:
            color = color_schema[allele]
        axes.bar(x_values, heights, bottom=bottom, width=width,
                 label=''.join(map(str, allele)), color=color)

        if bottom is None:
            bottom = heights
        else:
            bottom += heights

    matplotlib_support.set_x_ticks(x_values, sorted_pops, axes=axes, rotation=45)

    return {'alleles': alleles, 'x_values': x_values, 'sorted_pops': sorted_pops}


def get_chunk_with_n_uniq_haplos(variations, start_pos, num_desired_haplos):

    gts = variations[GT_FIELD]

    stop_pos = start_pos + 1

    n_haplos = 0
    last_n_vars = None
    while n_haplos < num_desired_haplos:
        this_gts = gts[start_pos:stop_pos, :, :]
        counts = _count_alleles(this_gts)
        n_haplos = len(counts)

        if n_haplos > num_desired_haplos:
            msg = 'Imposible to select a region with exactly the asked num of desired haplos'
            raise RuntimeError(msg)

        if n_haplos < num_desired_haplos:
            stop_pos += 1
            continue

        n_vars = this_gts.shape[0]
        if last_n_vars is None:
            last_n_vars = n_vars
        elif last_n_vars == n_vars:
            msg = 'There are no more vars left'
            raise RuntimeError(msg)
        else:
            last_n_vars = n_vars

        break

    return {'variations': variations.get_chunk(slice(start_pos, stop_pos)),
            'haplo_counts': counts}


def get_chunk_for_gene(variations, genes, gene_id, min_num_desired_haplos=0):
    gene = genes.get_gene(gene_id)

    vars_index = variations.pos_index
    gene = genes.get_gene(gene_id)
    idx0 = vars_index.index_pos(gene['Chromosome'].encode(), gene['Start'])
    idx1 = vars_index.index_pos(gene['Chromosome'].encode(), gene['End'])

    chunk = variations.get_chunk(slice(idx0, idx1 + 1))

    gts = chunk[GT_FIELD]
    counts = _count_alleles(gts)

    if len(counts) < min_num_desired_haplos:
        raise NotImplementedError('expand the region to get more haplos')

    return chunk


def _name_the_haplos(freqs, ref_pops):

    if freqs.shape[1] < len(ref_pops):
        msg = 'There are more ref pops than pops in freqs'
        raise ValueError(msg)

    most_freq_haplo_idx_for_pops = []
    for ref_pop in ref_pops:
        most_freq_haplo_idx_for_pop = numpy.argmax(freqs.loc[:, ref_pop].values)
        most_freq_haplo_idx_for_pops.append(most_freq_haplo_idx_for_pop)

    if len(most_freq_haplo_idx_for_pops) > len(set(most_freq_haplo_idx_for_pops)):
        msg = 'The most freq haplotype for some pops is the same one'
        raise RuntimeError(msg)

    new_haplo_names = []
    other_idx = 1
    for haplo_idx in range(freqs.shape[0]):
        if haplo_idx in most_freq_haplo_idx_for_pops:
            name = ref_pops[most_freq_haplo_idx_for_pops.index(haplo_idx)]
        else:
            name = f'other_{other_idx}'
            other_idx += 1
        new_haplo_names.append(name)
    freqs.index = new_haplo_names


def calc_mean_haplo_allele_freqs(variations, ref_pops, pops, n_succesful_attempts=100, n_max_attempts=None):

    if n_max_attempts is None:
        n_max_attempts = n_succesful_attempts * 100

    num_vars = variations.num_variations

    num_desired_haplos = len(ref_pops) + 1

    num_attempts = 0
    num_attemps_failed_due_to_n_haplos = 0
    num_attemps_failed_due_to_shared_haplos_between_ref_pops = 0
    freqss = []
    sorted_columns = None
    while len(freqss) < n_succesful_attempts:
        num_attempts += 1
        if num_attempts > n_max_attempts:
            raise RuntimeError('Num. max attempts exceeded')
        start_pos = random.randint(0, num_vars)

        try:
            res = get_chunk_with_n_uniq_haplos(variations, start_pos=start_pos, num_desired_haplos=num_desired_haplos)
        except RuntimeError:
            num_attemps_failed_due_to_n_haplos += 1
            continue

        freqs = calc_allele_freq(res['variations'], pops)

        try:
            _name_the_haplos(freqs, ref_pops=ref_pops)
        except RuntimeError:
            num_attemps_failed_due_to_shared_haplos_between_ref_pops += 1
            continue

        if sorted_columns is None:
            sorted_columns = sorted(freqs.columns, key=str)

        freqs = freqs.reindex(index=ref_pops + sorted(set(freqs.index).difference(ref_pops)),
                              columns=sorted_columns)

        freqss.append(freqs)

    total = None
    for freqs in freqss:
        if total is None:
            total = freqs
        else:
            total = total + freqs
    freqs = freqs / len(freqss)

    return {'mean_freqs': freqs, 'num_attemps_failed_due_to_n_haplos': num_attemps_failed_due_to_n_haplos,
            'num_attemps_failed_due_to_shared_haplos_between_ref_pops': num_attemps_failed_due_to_shared_haplos_between_ref_pops}


if __name__ == '__main__':

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, passports)

    res = get_chunk_with_n_uniq_haplos(variations, start_pos=80, num_desired_haplos=4)
    freqs = calc_allele_freq(res['variations'], pops)

    ref_pops = ['sp_pe', 'sp_ec' ,'sll_mx']
    _name_the_haplos(freqs, ref_pops=ref_pops)

    res = calc_mean_haplo_allele_freqs(variations, ref_pops, pops, n_succesful_attempts=100)
    freqs = res['mean_freqs']

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    plot_allele_freqs_bars(freqs, axes)

    plot_path = config.TMP_DIR / 'allele_freqs.svg'
    fig.tight_layout()
    fig.savefig(plot_path)
