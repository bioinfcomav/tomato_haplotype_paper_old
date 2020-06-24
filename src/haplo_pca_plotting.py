
from collections import defaultdict
from pprint import pprint

import numpy
import pandas

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import patches
from matplotlib import colors

from pca import write_curlywhirly_file
from colors import (ColorSchema, POP_COLORS, PINK_BLUE_CMAP_R, PINK_BLUE_CMAP_R2,
                    CLASSIFICATION_RANK1_COLORS, ELLIPSE_COLORS)
from haplo import parse_haplo_id
from haplo_pca import stack_aligned_pcas_projections



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


def write_pcas_curly_file(pcoas, out_dir, populations, haplo_classification=None):

    pop_for_samples = {sample: pop for pop, samples in populations.items() for sample in samples}

    path = out_dir / 'pcoas_along_the_genome.curly'

    all_projections = stack_aligned_pcas_projections(pcoas)
    pop_classification = {haplo_id: pop_for_samples[parse_haplo_id(haplo_id)[2]] for haplo_id in all_projections.index}

    categories = {'pop': pop_classification}

    if haplo_classification is not None:
        haplo_classification = {haplo_id: pop for pop, haplo_ids in haplo_classification.items() for haplo_id in haplo_ids}
        categories['haplo_class'] = haplo_classification

    write_curlywhirly_file(all_projections, path,
                           categories=categories)


def fit_ellipsoid(X, Y, scale):
    cov = numpy.cov(X, Y)
    eigvals, eigvecs = numpy.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    eigvecs = eigvecs.T

    width, height = 2 * numpy.sqrt(eigvals) * scale

    vx, vy = eigvecs[:,0][0], eigvecs[:,0][1]
    theta = numpy.arctan2(vy, vx)

    center = numpy.mean(X), numpy.mean(Y)

    ellipsoid = {'center': center, 'eigvals': eigvals, 'eigvecs': eigvecs, 'width': width, 'height': height, 'theta': theta}
    return ellipsoid


def plot_ellipsoids(axes, ellipsoids):
    if not ellipsoids:
        return

    for haplo_kind, ellipsoid in ellipsoids.items():
        color = ELLIPSE_COLORS.get(haplo_kind, '#d8d8d8')
        ellipse = patches.Ellipse(ellipsoid['center'],
                                  ellipsoid['width'],
                                  ellipsoid['height'],
                                  angle=numpy.degrees(ellipsoid['theta']),
                                  facecolor='none',
                                  edgecolor=color,
                                  linewidth=6,
                                  zorder=30,
                                  )
        axes.add_patch(ellipse)


def plot_hist2d_in_axes(aligned_pcoas_df, axes, x_lims=None, y_lims=None, ellipsoids=None):

    res = axes.hist2d(aligned_pcoas_df.values[:, 0], aligned_pcoas_df.values[:, 1], bins=100,
               norm=colors.LogNorm(), cmap=PINK_BLUE_CMAP_R2, zorder=10)

    plot_ellipsoids(axes, ellipsoids)

    if x_lims:
        axes.set_xlim(x_lims)
    if y_lims:
        axes.set_ylim(y_lims)
    return res


def plot_hist2d(aligned_pcoas_df, plot_path, x_lims=None, y_lims=None, ellipsoids=None):
    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    plot_hist2d_in_axes(aligned_pcoas_df, axes=axes,
                        x_lims=x_lims, y_lims=y_lims, ellipsoids=ellipsoids)

    fig.tight_layout()
    fig.savefig(str(plot_path))


def plot_classifications(aligned_pcoas_df, haplo_classification,
                         outlier_classes, plot_path, classes_to_ignore=None,
                         x_lims=None, y_lims=None, ellipsoids=None):

    if classes_to_ignore is None:
        classes_to_ignore = []

    fig = Figure((4, 15))
    FigureCanvas(fig) # Don't remove it or savefig will fail later

    all_haplo_ids_in_pca = set(aligned_pcoas_df.index)

    classes = sorted(haplo_classification.keys(), key=str)
    classes.sort(key=lambda x: x in outlier_classes)

    classes = [klass for klass in classes if klass not in classes_to_ignore]

    num_haplo_classes = len(classes)

    for klass_idx, klass in enumerate(classes):

        haplo_ids = haplo_classification[klass]

        haplo_ids = all_haplo_ids_in_pca.intersection(haplo_ids)

        axes = fig.add_subplot(num_haplo_classes, 1, klass_idx + 1)
        aligned_pcoas_df_for_this_class = aligned_pcoas_df.reindex(haplo_ids)
        x_values = aligned_pcoas_df_for_this_class.values[:, 0]
        y_values = aligned_pcoas_df_for_this_class.values[:, 1]
        axes.hist2d(x_values, y_values, bins=50,
                    norm=colors.LogNorm(), cmap=PINK_BLUE_CMAP_R2, zorder=10)
        plot_ellipsoids(axes, ellipsoids)

        if x_lims:
            axes.set_xlim(x_lims)
        if y_lims:
            axes.set_ylim(y_lims)

        axes.set_ylabel(klass)

    fig.tight_layout()
    fig.savefig(str(plot_path))


def calc_ellipsoids(classification, aligned_pcoas_df, scale=1, classes_to_ignore=None):
    classes = set(classification.values())

    if classes_to_ignore is not None:
        classes = classes.difference(classes_to_ignore)

    ellipsoids = {}
    for klass in classes:
        haplo_ids_in_class = [haplo_id for haplo_id in aligned_pcoas_df.index if classification[haplo_id] == klass]
        if not haplo_ids_in_class:
            continue
        haplo_positons_in_pcoas_for_class = aligned_pcoas_df.loc[haplo_ids_in_class, :]
        
        ellipsoid = fit_ellipsoid(haplo_positons_in_pcoas_for_class.iloc[:, 0],
                                  haplo_positons_in_pcoas_for_class.iloc[:, 1],
                                  scale=scale)
        ellipsoids[klass] = ellipsoid
    return ellipsoids


def filter_aligned_pcoas_df_for_samples(aligned_pcoas_df, samples):
    samples = set(samples)
    haplo_ids_to_keep = [haplo_id for haplo_id in aligned_pcoas_df.index if parse_haplo_id(haplo_id)[2] in samples]
    return aligned_pcoas_df.loc[haplo_ids_to_keep, :]


def _plot_pcas_for_samples(aligned_pcoas_df, samples, plot_path, color,
                           x_lims, y_lims, ellipsoids):

        aligned_pcoas_df_for_pop = filter_aligned_pcoas_df_for_samples(aligned_pcoas_df,
                                                                       samples)

        pop_fig = Figure((5, 10))
        FigureCanvas(pop_fig) # Don't remove it or savefig will fail later
        scatter_axes = pop_fig.add_subplot(211)

        hist2d_axes = pop_fig.add_subplot(212)

        x_values = aligned_pcoas_df_for_pop.values[:, 0]
        y_values = aligned_pcoas_df_for_pop.values[:, 1]

        scatter_axes.scatter(x_values, y_values, color=color, alpha=0.05, zorder=10)
        scatter_axes.set_facecolor('white')

        hist2d_axes.hist2d(x_values, y_values, bins=50,
                           norm=colors.LogNorm(), cmap=PINK_BLUE_CMAP_R2, zorder=10)

        plot_ellipsoids(scatter_axes, ellipsoids)
        plot_ellipsoids(hist2d_axes, ellipsoids)

        if x_lims:
            scatter_axes.set_xlim(x_lims)
            hist2d_axes.set_xlim(x_lims)
        if y_lims:
            scatter_axes.set_ylim(y_lims)
            hist2d_axes.set_ylim(y_lims)

        pop_fig.savefig(str(plot_path))
        scatter_axes.cla()
        hist2d_axes.cla()
        del scatter_axes
        del hist2d_axes
        pop_fig.clf()
        del pop_fig


def plot_pcas_per_pop(aligned_pcoas_df, out_dir,
                      populations, color_schema=None,
                      x_lims=None, y_lims=None,
                      ellipsoids=None):

    pop_names = sorted(populations.keys(), key=lambda x: '' if x is None else x)

    if color_schema is None:
        color_schema = ColorSchema(POP_COLORS)

    for idx, pop_name in enumerate(pop_names):

        samples_in_pop = populations[pop_name]

        color = color_schema[pop_name]
        pop_plot_path = out_dir / f'plots_for_pop_{pop_name}.png'

        _plot_pcas_for_samples(aligned_pcoas_df, samples=samples_in_pop,
                               plot_path=pop_plot_path,
                               color=color,
                               x_lims=x_lims, y_lims=y_lims,
                               ellipsoids=ellipsoids)

def plot_pcas_per_sample(aligned_pcoas_df, out_dir,
                         x_lims, y_lims,
                         populations,
                         ellipsoids=None,
                         color_schema=None):

    if color_schema is None:
        color_schema = ColorSchema(POP_COLORS)

    pops_for_samples = {sample: pop for pop, samples in populations.items() for sample in samples}

    fig_idx = 0
    for idx, (pop_name, samples) in enumerate(populations.items()):
        print(f'plotting procusted pcas for samples in pop: {pop_name}')

        color = color_schema[pop_name]

        pop_out_dir = out_dir / str(pop_name)
        pop_out_dir.mkdir(exist_ok=True)

        for sample in samples:
            plot_path = pop_out_dir / f'{sample}.svg'
            _plot_pcas_for_samples(aligned_pcoas_df, samples=[sample],
                                plot_path=plot_path,
                                color=color,
                                x_lims=x_lims, y_lims=y_lims,
                                ellipsoids=ellipsoids)
