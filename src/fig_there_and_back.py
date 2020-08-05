
import config

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5

import passport
import pop_building
from threre_and_back import calculate_and_plot_diversity_in_one_pop_vs_introgression_in_another
import matplotlib_support
import labels


if __name__ == '__main__':
    vars_path = config.WORKING_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops_rank1 = pop_building.get_pops(pops_descriptions, passports)
    pop_labels = labels.LABELS

    win_size = 10000

    founder_pop = 'slc_ma'

    samples_in_founder_pop = pops_rank1[founder_pop]

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    plot_path = config.FIG_THERE_AND_BACK
    diversity_index = 'num_poly95'

    introgression_source_pop = 'sp_pe'
    target_pop = 'slc_pe_n'
    label = pop_labels[introgression_source_pop] + ' → ' + pop_labels[target_pop]
    samples_in_pop_with_possible_introgressions = pops_rank1[target_pop]
    samples_in_introgression_source_pop = pops_rank1[introgression_source_pop]
    calculate_and_plot_diversity_in_one_pop_vs_introgression_in_another(variations,
                                                                        samples_in_pop_with_possible_introgressions,
                                                                        samples_in_founder_pop,
                                                                        samples_in_introgression_source_pop,
                                                                        diversity_index,
                                                                        win_size,
                                                                        axes,
                                                                        plot_line_label=label)

    introgression_source_pop = 'sp_pe'
    target_pop = 'slc_pe_s'
    label = pop_labels[introgression_source_pop] + ' → ' + pop_labels[target_pop]
    samples_in_pop_with_possible_introgressions = pops_rank1[target_pop]
    samples_in_introgression_source_pop = pops_rank1[introgression_source_pop]
    calculate_and_plot_diversity_in_one_pop_vs_introgression_in_another(variations,
                                                                        samples_in_pop_with_possible_introgressions,
                                                                        samples_in_founder_pop,
                                                                        samples_in_introgression_source_pop,
                                                                        diversity_index,
                                                                        win_size,
                                                                        axes,
                                                                        plot_line_label=label)

    introgression_source_pop = 'sp_ec'
    target_pop = 'slc_ec'
    label = pop_labels[introgression_source_pop] + ' → ' + pop_labels[target_pop]
    samples_in_pop_with_possible_introgressions = pops_rank1[target_pop]
    samples_in_introgression_source_pop = pops_rank1[introgression_source_pop]
    calculate_and_plot_diversity_in_one_pop_vs_introgression_in_another(variations,
                                                                        samples_in_pop_with_possible_introgressions,
                                                                        samples_in_founder_pop,
                                                                        samples_in_introgression_source_pop,
                                                                        diversity_index,
                                                                        win_size,
                                                                        axes,
                                                                        plot_line_label=label)


    matplotlib_support.set_axes_background(axes)
    axes.legend()
    axes.set_xlabel(f'Mean introgression freq. in Andean region')
    y_label = f'Mean num. poly. vars. (95%) ratio in {pop_labels[founder_pop]}'
    axes.set_ylabel(y_label)

    fig.tight_layout()
    fig.savefig(plot_path)

