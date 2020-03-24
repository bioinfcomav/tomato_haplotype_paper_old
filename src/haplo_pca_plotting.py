
from collections import defaultdict

import numpy
import pandas

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import patches

from pca import write_curlywhirly_file
from colors import (ColorSchema, POP_COLORS, PINK_BLUE_CMAP_R,
                    CLASSIFICATION_RANK1_COLORS)

ELLIPSE_COLORS = ColorSchema(CLASSIFICATION_RANK1_COLORS)


def plot_ellipsoids(axes, ellipsoids):
    if not ellipsoids:
        return

    for haplo_kind, ellipsoid in ellipsoids.items():
        color = ELLIPSE_COLORS[haplo_kind]
        ellipse = patches.Ellipse(ellipsoid['center'],
                                    ellipsoid['width'],
                                    ellipsoid['height'],
                                    angle=numpy.degrees(ellipsoid['theta']),
                                    facecolor='none',
                                    edgecolor=color,
                                    linewidth=2,
                                    zorder=1,
                                    alpha=0.7)
        axes.add_patch(ellipse)


def plot_haplo_pcas(pcoas, out_dir, populations=None, ellipsoids=None):

    if populations is None:
        there_are_pops = False
        populations = {None: sample for sample, _ in projections.index}
    else:
        there_are_pops = True

    pop_names = sorted(populations.keys(), key=lambda x: '' if x is None else x)

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    color_schema = ColorSchema(POP_COLORS)

    for pop_name in pop_names:
        samples_in_pop = populations[pop_name]

        for pcoa in pcoas:
            projections = pcoa['projections']

            xs = projections.iloc[:, 0]
            ys = projections.iloc[:, 1]

            if there_are_pops:
                haplo_is_in_pop = [sample in samples_in_pop for sample, _ in projections.index]
                xs_for_this_pop = xs[haplo_is_in_pop]
                ys_for_this_pop = ys[haplo_is_in_pop]

                color = color_schema[pop_name]
            else:
                color=None
                xs_for_this_pop = xs
                ys_for_this_pop = ys

            axes.scatter(xs_for_this_pop, ys_for_this_pop, color=color, alpha=0.05)

    x_lims =  axes.get_xlim()
    y_lims = axes.get_ylim()

    if ellipsoids:
        plot_ellipsoids(axes, ellipsoids)

    x_width = 1
    y_width = 1
    rectangles = []
    for idx, pop in enumerate(pop_names):
        color = color_schema[pop]
        rectangle = patches.Rectangle((idx, 0), x_width, y_width, color=color)
        rectangles.append(rectangle)
    axes.legend(rectangles, pop_names,
                fontsize='small')

    fig.tight_layout()

    plot_path = out_dir / 'haplo_pca.png'
    fig.savefig(str(plot_path))
    return {'x_lims': x_lims, 'y_lims': y_lims, 'color_schema': color_schema}


def plot_pcas_per_pop(pcoas, out_dir, x_lims, y_lims,
                      color_schema=None,
                      populations=None, ellipsoids=None):

    if populations is None:
        there_are_pops = False
        populations = {None: sample for sample, _ in projections.index}
    else:
        there_are_pops = True

    pop_names = sorted(populations.keys(), key=lambda x: '' if x is None else x)

    if color_schema is None:
        color_schema = ColorSchema(POP_COLORS)

    for idx, pop_name in enumerate(pop_names):

        samples_in_pop = populations[pop_name]

        pop_fig = Figure((5, 10))
        FigureCanvas(pop_fig) # Don't remove it or savefig will fail later
        pop_axes = pop_fig.add_subplot(211)
        pop_axes.set_xlim(x_lims)
        pop_axes.set_ylim(y_lims)

        hex_axes = pop_fig.add_subplot(212)
        hex_axes.set_xlim(x_lims)
        hex_axes.set_ylim(y_lims)

        for pcoa in pcoas:
            projections = pcoa['projections']

            if there_are_pops:
                haplo_is_in_pop = [sample in samples_in_pop for sample, _ in projections.index]
                xys_in_pop = projections.values[haplo_is_in_pop, :2]
            else:
                xys_in_pop = projections.values[:, :2]

            #seaborn.kdeplot(xys_in_pop[:, 0], xys_in_pop[:, 1], shade=False, cut=5, ax=pop_axes)
            color = color_schema[pop_name]
            pop_axes.scatter(xys_in_pop[:, 0], xys_in_pop[:, 1], color=color, alpha=0.05)
            pop_axes.set_facecolor('white')

            hex_axes.hexbin(xys_in_pop[:, 0], xys_in_pop[:, 1],
                            extent=(x_lims[0], x_lims[1], y_lims[0], y_lims[1]),
                            bins=None, gridsize=25, cmap=PINK_BLUE_CMAP_R)

        if ellipsoids:
            plot_ellipsoids(pop_axes, ellipsoids)

        pop_plot_path = out_dir / f'plots_for_pop_{pop_name}.png'

        pop_fig.savefig(str(pop_plot_path))
        pop_axes.cla()
        hex_axes.cla()
        del pop_axes
        del hex_axes
        pop_fig.clf()
        del pop_fig


def plot_pcas_per_sample(pcoas, out_dir,
                         x_lims, y_lims,
                         populations=None,
                         ellipsoids=None,
                         color_schema=None):

    if populations is None:
        there_are_pops = False
        populations = {None: sample for sample, _ in projections.index}
    else:
        there_are_pops = True

    pop_names = sorted(populations.keys(), key=lambda x: '' if x is None else x)

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    if color_schema is None:
        color_schema = ColorSchema(POP_COLORS)

    pops_for_samples = {sample: pop for pop, samples in populations.items() for sample in samples}

    samples_xs = {}
    samples_ys = {}
    for pcoa in pcoas:
        projections = pcoa['projections']
        xs = projections.iloc[:, 0]
        ys = projections.iloc[:, 1]
        for sample_idx, sample_haploid_idx in enumerate(projections.index):
            sample = sample_haploid_idx[0]
            pop = pops_for_samples[sample]
            if pop not in samples_xs:
                samples_xs[pop] = defaultdict(list)
                samples_ys[pop] = defaultdict(list)
            samples_xs[pop][sample].append(xs.iloc[sample_idx])
            samples_ys[pop][sample].append(ys.iloc[sample_idx])

    fig_idx = 0
    for pop in samples_xs.keys():
        print(f'plotting procusted pcas for samples in pop: {pop}')
        samples_xs_for_this_pop = samples_xs[pop]
        samples_ys_for_this_pop = samples_ys[pop]

        color = color_schema[pop]

        pop_out_dir = out_dir / str(pop)
        pop_out_dir.mkdir(exist_ok=True)

        for sample in samples_xs_for_this_pop.keys():

            fig = Figure((5, 10))
            fig_idx += 1
            FigureCanvas(fig) # Don't remove it or savefig will fail later
            scatter_axes = fig.add_subplot(211)
            scatter_axes.set_xlim(x_lims)
            scatter_axes.set_ylim(y_lims)
            scatter_axes.set_facecolor('white')

            hex_axes = fig.add_subplot(212)
            hex_axes.set_xlim(x_lims)
            hex_axes.set_ylim(y_lims)

            sample_xs = samples_xs_for_this_pop[sample]
            sample_ys = samples_ys_for_this_pop[sample]
            scatter_axes.scatter(sample_xs, sample_ys, color=color, alpha=0.1)

            if ellipsoids:
                plot_ellipsoids(scatter_axes, ellipsoids)

            hex_axes.hexbin(sample_xs, sample_ys,
                            extent=(x_lims[0], x_lims[1], y_lims[0], y_lims[1]),
                            bins=None, gridsize=25, cmap=PINK_BLUE_CMAP_R)
            
            plot_path = pop_out_dir / f'{sample}.svg'
            fig.savefig(str(plot_path))
            scatter_axes.cla()
            hex_axes.cla()
            del scatter_axes
            del hex_axes
            fig.clf()
            del fig
