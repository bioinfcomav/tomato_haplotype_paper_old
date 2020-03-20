
import config

from collections import defaultdict

import numpy
import pandas

from variation.variations import VariationsH5

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import patches

from passport import get_sample_passports
from haplo_pca import do_pcoas_along_the_genome
from procrustes import align_pcas_using_procrustes
from pca import write_curlywhirly_file
from pop_building import get_pops
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


def write_pcas_curly_file(pcoas, out_dir, populations=None, haplo_classification=None):

    pop_for_samples = {sample: pop for pop, samples in populations.items() for sample in samples}

    path = out_dir / 'pcoas_along_the_genome.curly'

    sample_names = []
    pops = []
    all_projections = []
    for pcoa_idx, pcoa in enumerate(pcoas):
        projections = pcoa['projections']
        chrom = pcoa['chrom'][-2:]

        if populations:
            sample_names.extend([f'{chrom}%{pcoa["win_start"]}%{sample}%{haploid_idx}' for sample, haploid_idx in projections.index])
            pops.extend([pop_for_samples[sample] for sample, _ in projections.index])

        all_projections.append(projections.values)

    all_projections = numpy.vstack(all_projections)

    categories = {}

    if populations:
        pop_classification = dict(zip(sample_names, pops))
        categories['pop'] = pop_classification

    if haplo_classification:
        haplo_classification = {haplo_id: pop for pop, haplo_ids in haplo_classification.items() for haplo_id in haplo_ids}
        categories['haplo_class'] = haplo_classification

    if not categories:
        categories = None

    all_projections = pandas.DataFrame(all_projections, index=sample_names)
    write_curlywhirly_file(all_projections, path,
                           categories=categories)


def get_sample_selection_criteria():

    rank1 = config.RANK1

    criteria = []
    samples_to_remove = []
    samples_to_keep = []

    return {'criteria': criteria, 'samples_to_remove': samples_to_remove,
            'samples_to_keep': samples_to_keep}


def get_haplotypes_to_exclude(path):
    haplos_to_exclude = defaultdict(list)

    if not path.exists():
        return {}

    for line in path.open('rt'):
        if line.startswith('Label'):
            continue
        line = line.strip()
        if not line:
            continue
        items = line.strip().split()[0].split('%')
        chrom = f'SL2.50ch{items[0]}'
        win_start = int(items[1])
        key = chrom, win_start
        sample = items[2]
        haploid_idx = int(items[3])
        haplos_to_exclude[key].append((sample, haploid_idx))
    return haplos_to_exclude


def get_haplotypes_to_include(path):
    return get_haplotypes_to_exclude(path)


if __name__ == '__main__':

    debug = True

    if debug:
        num_wins_to_process = 2
    else:
        num_wins_to_process = None

    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}

    vars_path = config.TIER1_PHASED_AND_IMPUTED_LOW_QUAL_09_MISSING_085
    variations = VariationsH5(str(vars_path), 'r')

    passports = get_sample_passports()

    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, passports)

    samples_to_use = {sample for samples in pops.values() for sample in samples}
    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        samples_to_use = samples_to_use.intersection(variations.samples)
    samples_to_use = sorted(samples_to_use)

    pcoas = do_pcoas_along_the_genome(variations, win_params,
                                      num_wins_to_process=num_wins_to_process,
                                      samples=samples_to_use, n_dims_to_keep=3)

    aligned_pcoas = list(align_pcas_using_procrustes(pcoas))

    out_dir = config.HAPLO_PCOA_DIR
    out_dir.mkdir(exist_ok=True)

    haplotypes_to_include_path = out_dir / 'haplotypes_to_include.txt'
    haplotypes_to_include = get_haplotypes_to_include(haplotypes_to_include_path)

    haplotypes_to_exclude_path = out_dir / 'haplotypes_to_exclude.txt'
    haplotypes_to_exclude = get_haplotypes_to_exclude(haplotypes_to_exclude_path)

    haplo_classification = None

    write_pcas_curly_file(aligned_pcoas, out_dir, pops, haplo_classification)

    ellipsoids = None

    res = plot_haplo_pcas(aligned_pcoas, out_dir, populations=pops,
                          ellipsoids=ellipsoids)

    per_pop_out_dir = out_dir / 'per_pop'
    per_pop_out_dir.mkdir(exist_ok=True)
    plot_pcas_per_pop(aligned_pcoas, per_pop_out_dir, res['x_lims'], res['y_lims'],
                      color_schema=res['color_schema'],
                      populations=pops, ellipsoids=ellipsoids)

    per_sample_out_dir = out_dir / 'per_sample'
    per_sample_out_dir.mkdir(exist_ok=True)
    plot_pcas_per_sample(aligned_pcoas, per_sample_out_dir,
                         res['x_lims'], res['y_lims'],
                         populations=pops,
                         ellipsoids=ellipsoids)
