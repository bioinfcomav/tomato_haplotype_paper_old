
import config

import math

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5

import passport
import pop_building
import matplotlib_support
import labels
from fig_haplo_pca import plot_haplo_pca_for_sample, HAPLO_PCOAS_X_LIMS, HAPLO_PCOAS_Y_LIMS


def plot_haplo_pcas_per_sample(fig, imputed_variations, samples_to_use, pops, cache_dir):

    pops_by_sample = {sample: pop for pop, samples in pops.items() for sample in samples}

    ncols = 10
    nrows = math.ceil(len(samples_to_use) / ncols)

    axes_col_widths = [1 / ncols] * ncols
    axes_row_heights = [1 / nrows] * nrows

    max_num_samples = None

    for sample_idx, sample in enumerate(sorted(samples_to_use)):
        print(sample)
        row_idx = sample_idx // ncols
        col_idx = sample_idx % ncols
        sample_axes = matplotlib_support.add_axes(fig, row_idx, col_idx,
                                                  axes_col_widths=axes_col_widths,
                                                  axes_row_heights=axes_row_heights,
                                                  top_margin=0.15)
        plot_haplo_pca_for_sample(sample_axes, imputed_variations, sample,
                                  samples_to_use, pops, cache_dir, marker='.')
        sample_axes.set_xlim(HAPLO_PCOAS_X_LIMS)
        sample_axes.set_ylim(HAPLO_PCOAS_Y_LIMS)

        matplotlib_support.set_axes_background(sample_axes)

        pop = labels.LABELS[pops_by_sample[sample]]
        sample_axes.set_title(f'{sample} ({pop})', fontdict={'fontsize':5})

        if row_idx != nrows:
            matplotlib_support.turn_off_x_axis(sample_axes)
        if col_idx != 0:
            matplotlib_support.turn_off_y_axis(sample_axes)
        if max_num_samples and sample_idx >= max_num_samples:
            break


if __name__ == '__main__':

    cache_dir = config.CACHE_DIR
    dist_threshold = config.CLASSIFICATION_CONFIG['thinning_dist_threshold']

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    imputed_variations = VariationsH5(str(vars_path), 'r')

    plot_path = config.FIG_HAPLO_PCA_PER_SAMPLE

    sample_passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, sample_passports)
    samples_to_use = {sample for samples in pops.values() for sample in samples}

    fig = Figure((20, 70))
    FigureCanvas(fig) # Don't remove it or savefig will fail later

    plot_haplo_pcas_per_sample(fig, imputed_variations, samples_to_use, pops, cache_dir)

    fig.savefig(plot_path)
