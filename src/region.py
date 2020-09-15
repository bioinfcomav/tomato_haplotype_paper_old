
import config

import math

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation import POS_FIELD
from variation.variations.filters import SNPPositionFilter, FLT_VARS

import imputation
import allele_freq
import haplo_net
import colors
import matplotlib_support
import rarefaction_vars


def plot_snp_density(axes, start, end, poss, n_bins=None):

    desired_mean_n_vars_per_bin = 5

    if n_bins is None:
        n_vars = poss.shape[0]
        n_bins = int(n_vars / desired_mean_n_vars_per_bin)
        if n_bins > 40:
            n_bins = 40

    bin_edges = numpy.linspace(start, end, n_bins + 1, dtype=int)

    n_vars_per_bin = []
    for bin_start, bin_end in zip(bin_edges[:-1], bin_edges[1:]):
        mask = numpy.logical_and(poss > bin_start, poss <= bin_end)
        n_vars_in_bin = numpy.sum(mask)
        n_vars_per_bin.append(n_vars_in_bin)

    x_poss = (bin_edges[:-1] + bin_edges[1:]) / 2
    width = x_poss[1] - x_poss[0]
    axes.bar(x_poss, n_vars_per_bin, width=width)
    axes.set_xlim((start, end))
    axes.set_ylabel('Num. vars.')
    return {'bin_edges': bin_edges}


def plot_mean_quantitative_index(axes, bin_edges, poss, values):

    means = []
    for bin_start, bin_end in zip(bin_edges[:-1], bin_edges[1:]):
        mask = numpy.logical_and(poss > bin_start, poss <= bin_end)
        this_bin_values = values[mask]
        if this_bin_values.size:
            mean = numpy.nanmean(this_bin_values)
        else:
            mean = 0
        means.append(mean)
    x_poss = (bin_edges[:-1] + bin_edges[1:]) / 2
    axes.scatter(x_poss, means)


def _plot_gene_arrows(axes, genes, y_pos):
    for gene in genes:
        axes.arrow(gene['Start'], y_pos, gene['End'] - gene['Start'], 0, width=0.1)
        axes.text(gene['Start'], y_pos - 0.1, gene['name'])


def plot_genes(axes, genes):

    fwd_genes = []
    rev_genes = []
    for gene in genes:
        if gene['Start'] < gene['End']:
            fwd_genes.append(gene)
        else:
            rev_genes.append(gene)

    _plot_gene_arrows(axes, fwd_genes, 1)
    _plot_gene_arrows(axes, rev_genes, 0.5)
    axes.set_ylim((0.3, 1.2))
    matplotlib_support.set_axes_background(axes)
    matplotlib_support.turn_off_y_axis(axes)


def analyze_region(variations, region_to_analyze, pops, pop_order, out_dir, genes, min_num_vars=0,
                   min_num_var_for_net=14, remove_indels=False, ignore_haplonet_failure=False):

    chrom = region_to_analyze['chrom']
    start = region_to_analyze['start']
    end = region_to_analyze['end']
    variations_in_region = get_variations_in_region(variations,
                                                    chrom,
                                                    start, end,
                                                    min_num_vars=min_num_vars)

    allowed_missing_gts = len(variations_in_region.samples) - 10
    res = rarefaction_vars.calc_pop_stats_per_var(variations_in_region, allowed_missing_gts=allowed_missing_gts)

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(311)
    res2 = plot_snp_density(axes, start, end,
                            variations_in_region[POS_FIELD])
    bin_edges = res2['bin_edges']

    axes = fig.add_subplot(312, sharex=axes)
    plot_mean_quantitative_index(axes, bin_edges, variations_in_region[POS_FIELD],
                                 res['unbiased_exp_het'])
    axes.set_ylabel('Exp. het.')

    axes = fig.add_subplot(313, sharex=axes)
    plot_genes(axes, genes.get_genes_in_region(chrom, start, end))

    plot_path = out_dir / f'snp_density.svg'
    fig.tight_layout()
    fig.savefig(plot_path)

    variations = imputation.impute_variations(variations_in_region,
                                              genetic_map_path=config.BEAGLE_MAP)

    freqs = allele_freq.calc_allele_freq(variations, pops)

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    allele_freq.plot_allele_freqs_bars(freqs, axes, pop_order=pop_order)

    plot_path = out_dir / f'allele_freqs.svg'
    fig.tight_layout()
    fig.savefig(plot_path)

    haplo_net_failed = False
    try:
        res = haplo_net.calc_haplotype_network(variations,
                                            min_num_snp_for_network=min_num_var_for_net,
                                            remove_indels=remove_indels)
    except:
        haplo_net_failed = True

    if not ignore_haplonet_failure and haplo_net_failed:
        raise RuntimeError('Haplonet failed')

    if not haplo_net_failed:
        pop_for_samples = {sample: pop for pop, samples in pops.items() for sample in samples}
        fig = Figure()
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)
        haplo_net.plot_network(res['network_nx'], res['sample_sets'],
                            axes, pop_for_samples,
                            pop_color_schema=colors.ColorSchema(colors.CLASSIFICATION_RANK1_COLORS))

        plot_path = out_dir / 'network.svg'
        fig.tight_layout()
        fig.savefig(plot_path)


def get_variations_in_region(variations, chrom, start, end, min_num_vars=0):

    chrom = chrom.encode()
    index = variations.pos_index
    idx0 = index.index_pos(chrom, start)
    idx1 = index.index_pos(chrom, end)

    extended = False
    if min_num_vars:
        num_vars = idx1 - idx0
        if num_vars < min_num_vars:
            vars_to_extend = math.ceil((min_num_vars - num_vars) / 2)
            idx0 = idx0 - vars_to_extend
            idx1 = idx1 + vars_to_extend
            extended = True
        
    chunk = variations.get_chunk(slice(idx0, idx1))

    if extended:
        chunk = SNPPositionFilter([(chrom, )])(chunk)[FLT_VARS]

    return chunk
