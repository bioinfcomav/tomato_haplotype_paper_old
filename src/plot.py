
import math

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import seaborn

from variation import CHROM_FIELD, POS_FIELD

from genome_coord_transform import (get_genome_sizes, GenomeCoordinateConverter,
                                    PositionInPericentromericRegion,
                                    GenomeCoordinateConverter2)
import colors


def get_sorted_legend_handles(axes):

    handles, labels = axes.get_legend_handles_labels()
    # sort both labels and handles by labels
    handles, labels = zip(*sorted(zip(handles, labels), key=lambda t: t[1]))
    return handles, labels


def plot_var_density_along_genome(variations, axes=None, plot_path=None,
                                  coord_converter=None,
                                  window_size=1e6, marker='.', linestyle='None',
                                  alpha=1):
    if axes is None and plot_path is None:
        raise ValueError('An axes or a plot path should be provided')
    if axes is not None and plot_path is not None:
        raise ValueError('Either axes or a plot path should be provided')

    if axes is None:
        fig = Figure((10, 5))
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)

    if coord_converter is None:
        chrom_lens = get_genome_sizes()
        coord_converter = GenomeCoordinateConverter(chrom_lens)

    converted_poss = []
    for chrom, pos in zip(variations[CHROM_FIELD], variations[POS_FIELD]):
        try:
            pos = coord_converter.transform_coordinate(chrom, pos)
        except PositionInPericentromericRegion:
            continue
        converted_poss.append(pos)
    converted_poss = numpy.array(converted_poss)
    
    min_ = min(converted_poss)
    max_ = max(converted_poss)
    n_bins = math.ceil((max_ - min_) / window_size)

    counts, bin_edges = numpy.histogram(converted_poss, bins=n_bins)

    xs = (bin_edges[:-1] + bin_edges[1:]) / 2

    poss_with_counts = counts != 0
    counts = counts[poss_with_counts]
    xs = xs[poss_with_counts]

    axes.plot(xs, counts, marker=marker, linestyle=linestyle)

    axes.set_ylim(0, axes.get_ylim()[1])

    draw_chromosome_lines(axes, coord_converter)

    if plot_path:
        fig.tight_layout()
        fig.savefig(str(plot_path))


def draw_chromosome_lines(axes, coord_converter, min_y=None, max_y=None):

    y_lims = axes.get_ylim()

    if min_y is None:
        min_y = y_lims[0]
    if max_y is None:
        max_y = y_lims[1]

    try:
        chrom_lens = coord_converter.chrom_lens
    except AttributeError:
        chrom_lens = None

    if chrom_lens:
        for chrom, length in chrom_lens.items():
            chrom_end = coord_converter.transform_coordinate(chrom, 1) + length
            axes.plot([chrom_end, chrom_end], [min_y, max_y],
                      c=colors.BLACK, linestyle='solid')

    try:
        pericentromeric_starts = coord_converter.pericentromeric_starts
    except AttributeError:
        pericentromeric_starts = None

    if pericentromeric_starts:
        for chrom, position in pericentromeric_starts.items():
            transformed_position = coord_converter.transform_coordinate(chrom, position)
            axes.plot([transformed_position, transformed_position],
                       [min_y, max_y], c=colors.DARK_GRAY)

    xticks = []
    xtick_labels = []
    for chrom, span in coord_converter.chrom_spans.items():
        xticks.append(span[1])
        xtick_labels.append(chrom)
    axes.set_xticklabels(xtick_labels, rotation=45, ha='right')
    axes.set_xticks(xticks)

    axes.grid(False)


def plot_histogram(vector, plot_path, labels=None, kde=False, line=True, axis_lims=None, n_bins=20):

    if labels is None:
        labels = {}

    vector = numpy.array(vector)
    vector = vector[numpy.logical_not(numpy.isnan(vector))]

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    range_ = None
    hist_kws = {}
    if axis_lims:
        if 'x_lim' in axis_lims:
            range_ = axis_lims['x_lim']
            hist_kws = {"range": range_}

    seaborn.distplot(vector, ax=axes, kde=kde, bins=n_bins, hist_kws=hist_kws)

    if line:
        if range_:
            counts, bin_edges = numpy.histogram(vector, bins=n_bins, range=range_)
        else:
            counts, bin_edges = numpy.histogram(vector, bins=n_bins)
        x_values = (bin_edges[:-1] + bin_edges[1:]) / 2
        axes.plot(x_values, counts)

    if 'x_axis' in labels:
        axes.set_xlabel(labels['x_axis'])
    if 'y_axis' in labels:
        axes.set_ylabel(labels['y_axis'])

    if axis_lims:
        if 'x_lim' in axis_lims:
            axes.set_xlim(axis_lims['x_lim'])
        if 'y_lim' in axis_lims:
            axes.set_ylim(axis_lims['y_lim'])

    fig.tight_layout()
    fig.savefig(str(plot_path))


def plot_scatter(x_values, y_values, plot_path, labels=None, fit_reg=False,
                 point_labels=None, color=None, alpha=1):

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    axes.set_facecolor('w')
    axes.grid(False)
    axes.spines['left'].set_linestyle("-")
    axes.spines['left'].set_color("black")
    axes.spines['bottom'].set_linestyle("-")
    axes.spines['bottom'].set_color("black")

    if color:
        seaborn.regplot(x_values, y_values, ax=axes, fit_reg=fit_reg, color=color, scatter_kws={'alpha':alpha})
    else:
        seaborn.regplot(x_values, y_values, ax=axes, fit_reg=fit_reg, scatter_kws={'alpha':alpha})

    if point_labels is not None:
        for label, x, y in zip(point_labels, x_values, y_values):
            axes.annotate(label, (x, y))

    if labels:
        if 'x_axis' in labels:
            axes.set_xlabel(labels['x_axis'])
        if 'y_axis' in labels:
            axes.set_ylabel(labels['y_axis'])

    axes.set_ylim((numpy.amin(y_values), numpy.amax(y_values)))
    axes.set_xlim((numpy.amin(x_values), numpy.amax(x_values)))

    fig.tight_layout()
    fig.savefig(str(plot_path))
    return {'axes': axes, 'fig': fig}
