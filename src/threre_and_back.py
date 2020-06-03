
import config

from array import array

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation import CHROM_FIELD, POS_FIELD
from variation.variations import VariationsH5
from variation.variations.filters import SampleFilter, FLT_VARS

import passport
import pop_building
import introgressions
from snp_filtering import keep_variations_variable_in_samples
from rarefaction_vars import calc_pop_stats_per_var
import matplotlib_support


N_BINS = 6


def calc_mean_per_window(values, chroms, poss, win_size, mean_operation, min_num_values=None):
    chrom_names = sorted(numpy.unique(chroms))
    chrom_sizes = {chrom: numpy.max(poss[chroms == chrom]) for chrom in chrom_names}

    means = array('f')
    for chrom, chrom_size in chrom_sizes.items():
        bins = numpy.arange(0, chrom_size, win_size)
        for start, stop in zip(bins[:-1], bins[1:]):
            in_bin = numpy.logical_and(chroms == chrom,
                                       numpy.logical_and(poss >= start, poss < stop))
            if min_num_values and numpy.sum(in_bin) < min_num_values:
                continue
            
            values_in_bin = values[in_bin]
            if min_num_values and values_in_bin.size < min_num_values:
                continue
            if mean_operation == 'mean':
                mean = numpy.nanmean(values_in_bin)
            elif mean_operation == 'ratio':
                if values_in_bin.size:
                    mean = numpy.sum(values_in_bin)  / values_in_bin.size
                else:
                    mean = numpy.nan
            means.append(mean)
    means = numpy.array(means)
    return {'means': means}


def _remove_items_with_nans(vector1, vector2):
    mask = numpy.logical_and(numpy.logical_not(numpy.isnan(vector1)),
                             numpy.logical_not(numpy.isnan(vector2)))
    return vector1[mask], vector2[mask]


def _plot_diversity_vs_introgression_freq(diversity_in_founder_pop, introgression_freq,
                                          axes, axes_num_vars=None, plot_type=None,
                                          plot_line_label=None,
                                          n_bins=10, introgression_freq_range=None,
                                          flier_symbol='', whisker_quartiles=(15, 85)):

    diversity_in_founder_pop, introgression_freq = _remove_items_with_nans(diversity_in_founder_pop, introgression_freq)

    if introgression_freq_range is None:
        introgression_freq_range = (numpy.min(introgression_freq), numpy.max(introgression_freq))

    bin_edges = numpy.linspace(introgression_freq_range[0], introgression_freq_range[1], num=n_bins + 1)
    x_poss = (bin_edges[:-1] + bin_edges[1:]) / 2
    width = (bin_edges[1] - bin_edges[0]) * 0.9

    means = []
    for min_intro_freq, max_intro_freq, x_pos in zip(bin_edges[:-1], bin_edges[1:], x_poss):
        mask = numpy.logical_and(introgression_freq > min_intro_freq,
                                 introgression_freq <= max_intro_freq)
        diversities_in_bin = diversity_in_founder_pop[mask]
        num_vars = numpy.sum(mask)

        if plot_type == 'boxplot':
            axes.boxplot(diversities_in_bin, positions=[x_pos], widths=[width],
                         sym=flier_symbol, whis=whisker_quartiles)
        if plot_type == 'mean':
            means.append(numpy.mean(diversities_in_bin))
        elif plot_type == 'bar':
            ratio_of_true_values = numpy.sum(diversities_in_bin) / diversities_in_bin.size
            axes.bar([x_pos], ratio_of_true_values, width=width)
        if axes_num_vars:
            axes_num_vars.bar([x_pos],num_vars, width=width)
    if plot_type == 'mean':
        axes.plot(x_poss, means, label=plot_line_label)
        

def plot_diversity_box_plot_vs_introgression_freq(diversity_in_founder_pop, introgression_freq,
                                                  axes, axes_num_vars,
                                                  n_bins=N_BINS, introgression_freq_range=None):
    _plot_diversity_vs_introgression_freq(diversity_in_founder_pop, introgression_freq,
                                          axes, axes_num_vars, plot_type='boxplot',
                                          n_bins=n_bins, introgression_freq_range=introgression_freq_range)


def plot_diversity_mean_vs_introgression_freq(diversity_in_founder_pop, introgression_freq,
                                                  axes, axes_num_vars, plot_line_label=None,
                                                  n_bins=N_BINS, introgression_freq_range=None):
    _plot_diversity_vs_introgression_freq(diversity_in_founder_pop, introgression_freq,
                                          axes, axes_num_vars, 
                                          plot_line_label=plot_line_label,
                                          plot_type='mean',
                                          n_bins=n_bins, introgression_freq_range=introgression_freq_range)


def plot_num_poly_vs_introgression_freq(var_is_poly95_in_founder_pop, introgression_freq,
                                        axes, axes_num_vars, plot_line_label=None,
                                        n_bins=N_BINS, introgression_freq_range=None):

    _plot_diversity_vs_introgression_freq(var_is_poly95_in_founder_pop, introgression_freq,
                                          axes, axes_num_vars,
                                          plot_line_label=plot_line_label,
                                          plot_type='bar',
                                          n_bins=n_bins, introgression_freq_range=introgression_freq_range)


def calculate_and_plot_diversity_in_one_pop_vs_introgression_in_another(variations,
                                                                        samples_in_pop_with_possible_introgressions,
                                                                        samples_in_founder_pop,
                                                                        samples_in_introgression_source_pop,
                                                                        diversity_index,
                                                                        win_size,
                                                                        axes,
                                                                        axes2=None,
                                                                        plot_line_label=None):

    all_samples = samples_in_pop_with_possible_introgressions + samples_in_founder_pop + samples_in_introgression_source_pop
    variations = keep_variations_variable_in_samples(variations, all_samples)

    vars_in_founder_pop = SampleFilter(samples_in_founder_pop)(variations)[FLT_VARS]

    allowed_missing_gts = len(samples_in_founder_pop) - 10
    res = calc_pop_stats_per_var(vars_in_founder_pop, allowed_missing_gts=allowed_missing_gts,
                                 stats_to_calc=['unbiased_exp_het', 'var_is_poly95'])
    nei_diversity_in_founder_pop = res['unbiased_exp_het']
    var_is_poly95_in_founder_pop = res['var_is_poly95']


    chroms = variations[CHROM_FIELD]
    poss = variations[POS_FIELD]
    if diversity_index == 'num_poly95':
        mean_diversity_values_per_win = calc_mean_per_window(var_is_poly95_in_founder_pop, chroms, poss, win_size, mean_operation='mean')['means']
    elif diversity_index == 'unbiased_exp_het':
        mean_diversity_values_per_win = calc_mean_per_window(nei_diversity_in_founder_pop, chroms, poss, win_size, mean_operation='ratio')['means']

    introgession_config = {'samples_in_pop_with_introgressions': samples_in_pop_with_possible_introgressions,
                           'samples_in_founder_pop': samples_in_founder_pop,
                           'samples_in_introgression_source_pop': samples_in_introgression_source_pop,
                           'freq_threshold_to_consider_allele_present_in_founder_pop' : config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_PRESENT,
                           'freq_threshold_to_consider_allele_common_in_introgression_source_pop': config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_COMMON,
                          }
    introgression_freq = introgressions.calc_introgession_freq_for_vars(variations, introgession_config)
    mean_introgression_freq_per_win = calc_mean_per_window(introgression_freq, chroms, poss, win_size, mean_operation='ratio')['means']

    plot_diversity_mean_vs_introgression_freq(mean_diversity_values_per_win,
                                                  mean_introgression_freq_per_win,
                                                  axes, axes2, plot_line_label=plot_line_label)

    if diversity_index == 'num_poly95':
        y_label = f'Mean num. poly. vars. (95%) ratio in win ({win_size} bp)'
    elif diversity_index == 'unbiased_exp_het':
        y_label = f'Mean unbiased expected het. in win ({win_size} bp)'

    axes.set_ylabel(y_label)
    if axes2:
        axes2.set_ylabel('Num. windows')
        matplotlib_support.turn_off_x_axis(axes)
        axes2.set_xlabel(f'Mean introgression freq. in win ({win_size} bp)')
    else:
        axes.set_xlabel(f'Mean introgression freq. in win ({win_size} bp)')

    #loc = plticker.MultipleLocator(base=0.2) # this locator puts ticks at regular intervals
    #loc = plticker.AutoLocator()
    #axes2.xaxis.set_major_locator(loc)
    

if __name__ == '__main__':
    vars_path = config.WORKING_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops_rank1 = pop_building.get_pops(pops_descriptions, passports)

    pops_descriptions = {config.RANK2: config.ALL_POPS}
    pops_rank2 = pop_building.get_pops(pops_descriptions, passports)

    win_size = 10000

    founder_pop = 'slc_ma'

    target_and_introgression_source = [('slc_pe', config.RANK1, 'sp_pe'),
                                       ('slc_ec', config.RANK1, 'sp_ec'),
                                       #('slc_pe_n_400m', config.RANK2, 'sp_pe'),
                                       #('slc_pe_n_1000m', config.RANK2, 'sp_pe'),
                                       #('slc_ec_n_600m', config.RANK2, 'sp_ec'),
                                       #('slc_ec_c_800m', config.RANK2, 'sp_ec')
                                       ]
    for target_pop, target_pop_rank, introgression_source_pop in target_and_introgression_source:
        if target_pop_rank == config.RANK1:
            samples_in_pop_with_possible_introgressions = pops_rank1[target_pop]
        else:
            samples_in_pop_with_possible_introgressions = pops_rank2[target_pop]

        samples_in_founder_pop = pops_rank1[founder_pop]
        samples_in_introgression_source_pop = pops_rank1[introgression_source_pop]

        out_dir = config.THERE_AND_BACK_DIR
        out_dir.mkdir(exist_ok=True)

        if True:
            plot_path = out_dir / f'num_poly95_in_{founder_pop}_vs_introgressions_in_{target_pop}_coming_from_{introgression_source_pop}.svg'
            diversity_index = 'num_poly95'
        else:
            plot_path = out_dir / f'unbiased_exp_het_in_{founder_pop}_vs_introgressions_in_{target_pop}_coming_from_{introgression_source_pop}.svg'
            diversity_index = 'unbiased_exp'

        fig = Figure((7, 10))
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(211)
        axes2 = fig.add_subplot(212, sharex=axes)
        axes2.set_yscale('log')

        calculate_and_plot_diversity_in_one_pop_vs_introgression_in_another(variations,
                                                                            samples_in_pop_with_possible_introgressions,
                                                                            samples_in_founder_pop,
                                                                            samples_in_introgression_source_pop,
                                                                            diversity_index,
                                                                            win_size,
                                                                            axes, axes2)
        matplotlib_support.set_axes_background(axes)
        matplotlib_support.set_axes_background(axes2)
        fig.tight_layout()
        fig.savefig(plot_path)
