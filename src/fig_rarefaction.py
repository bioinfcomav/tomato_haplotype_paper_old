
import config

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5

import passport
import pop_building
import colors
import labels
import rarefaction_vars
import matplotlib_support


def _plot_rarefacted_diversities_for_pops(diversities_per_pop, num_samples_per_pop, axes, pop_colors=None, y_lims=None):
    pops = sorted(diversities_per_pop.keys(), key=str)

    color_schema = colors.ColorSchema(pop_colors)

    for pop in pops:
        color = color_schema[pop]
        x_values = num_samples_per_pop[pop]
        pop_values_for_diversity_field = diversities_per_pop[pop]

        if isinstance(pop_values_for_diversity_field[0], (float, int)):
            y_values = pop_values_for_diversity_field
            axes.plot(x_values, y_values, label=labels.LABELS[pop], color=color)
        elif isinstance(pop_values_for_diversity_field[0], numpy.ndarray) and pop_values_for_diversity_field[0].size == 3:
            quartiles_low = [quartile_values[0] for quartile_values in pop_values_for_diversity_field]
            medians = [quartile_values[1] for quartile_values in pop_values_for_diversity_field]
            quartiles_high = [quartile_values[-1] for quartile_values in pop_values_for_diversity_field]
            axes.plot(x_values, medians, label=pop, color=color, zorder=10)
            axes.fill_between(x_values, quartiles_low, quartiles_high, alpha=0.5, color=color, zorder=5)

    if y_lims is not None:
        axes.set_ylim(*y_lims)


def plot_rarefacted_diversities(rarefacted_diversities, fig, pop_sets,
                                diversity_indexes=None, pop_colors=None, y_lims=None):

    if diversity_indexes is None:
        diversity_indexes = sorted({diversity_index for diversities in rarefacted_diversities.values() for diversity_index in diversities.keys()})

    num_diversity_indexes = len(diversity_indexes)
    num_pop_sets = len(pop_sets)

    axes_col_widths = [1 / num_pop_sets] * num_pop_sets
    axes_row_heights = [1 / num_diversity_indexes] * num_diversity_indexes

    for pop_set_idx, pop_set in enumerate(pop_sets):
        for diversity_index_idx, diversity_index in enumerate(diversity_indexes):
            col_idx = pop_set_idx
            row_idx = diversity_index_idx
            axes = matplotlib_support.add_axes(fig, row_idx=row_idx,
                                               col_idx=col_idx,
                                               axes_row_heights=axes_row_heights, axes_col_widths=axes_col_widths)
            diversities_per_pop = {pop: rarefacted_diversities[pop][diversity_index] for pop in pop_set}
            num_samples_per_pop = {pop: rarefacted_diversities[pop]['num_samples'] for pop in pops}
        
            _plot_rarefacted_diversities_for_pops(diversities_per_pop, num_samples_per_pop, axes, pop_colors=pop_colors, y_lims=y_lims.get(diversity_index))

            matplotlib_support.set_axes_background(axes)

            if not col_idx == 0:
                matplotlib_support.turn_off_y_axis(axes)
            else:
                axes.set_ylabel(diversity_index)
            if row_idx != (num_diversity_indexes - 1):
                matplotlib_support.turn_off_x_axis(axes)
            else:
                axes.set_xlabel('Num. accs.')
                axes.legend(loc='upper right')


if __name__ == '__main__':

    rarefaction_range = (8, 23)
    #rarefaction_range = (8, 10)

    vars_path = config.WORKING_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()

    main_pops = ['sp_pe' ,'sp_ec', 'slc_ma', 'slc_ec', 'slc_pe', 'sll_mx']
    vintage_pops = ['sll_vint', 'slc_world']
    hybrid_pops = ['sll_modern', 'sp_x_sl', 'sp_x_sp']
    all_pops = main_pops + vintage_pops + hybrid_pops
    pop_sets = [main_pops, vintage_pops, hybrid_pops]

    pops_descriptions = {config.RANK1: all_pops}
    pops = pop_building.get_pops(pops_descriptions, passports)
    pop_colors = colors.CLASSIFICATION_RANK1_COLORS

    rarefacted_diversities = rarefaction_vars.calc_rarefacted_diversities(variations, pops,
                                                                          rarefaction_range)

    y_lims = {'mean_num_alleles': (1, 1.6),
              'poly80': (0, 20),
              'poly95': (0, 60),
              'ratio_poly80/poly95': (0, 0.6),
              'unbiased_exp_het_percentiles': (0, 0.35),
              'unbiased_exp_het_mean': (0, 0.20)
              }
    diversity_indexes = ['unbiased_exp_het_mean', 'poly95', 'mean_num_alleles']
    #diversity_indexes=None

    fig = Figure((10, 7))
    FigureCanvas(fig) # Don't remove it or savefig will fail later

    plot_rarefacted_diversities(rarefacted_diversities, fig, pop_colors=pop_colors, pop_sets=pop_sets, diversity_indexes=diversity_indexes, y_lims=y_lims)

    fig.savefig(config.FIG_RAREFACTION)
