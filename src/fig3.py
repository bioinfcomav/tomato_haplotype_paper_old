
import config

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import cartopy.crs as ccrs

from variation.variations import VariationsH5

from haplo_auto_classification import (detected_outliers_and_classify_haplos,
                                       calc_haplo_pop_composition_freq_dframe)
import passport
from pop_building import get_pops
from haplo_pca_plotting import filter_aligned_pcoas_df_for_samples
import colors
import labels
import matplotlib_support
from fig1 import HAPLO_PCOAS_X_LIMS, HAPLO_PCOAS_Y_LIMS
from plot_geo_map import plot_geo_rank1_for_main_pops
from faststructure_parse_results import parse_faststructure_results
from faststructure_plot import _calc_admixtures_per_pop

def plot_haplo_pca_for_samples(aligned_pcoas_df, haplo_classification,
                               samples, axes,
                               color_schema=None,
                               alpha=1):
    aligned_pcoas_df_for_sample = filter_aligned_pcoas_df_for_samples(aligned_pcoas_df,
                                                                      samples)
    
    haplo_classes = numpy.array([haplo_classification[haplo_id] for haplo_id in aligned_pcoas_df_for_sample.index])
    uniq_haplo_classes = sorted(set(haplo_classes))
    
    for haplo_class in uniq_haplo_classes:
        mask = haplo_classes == haplo_class
        aligned_pcoas_for_sample_and_haplo_class = aligned_pcoas_df_for_sample.values[mask, :]
        color = colors.HAPLO_COLORS.get(haplo_class, None)
        axes.scatter(aligned_pcoas_for_sample_and_haplo_class[:, 0],
                     aligned_pcoas_for_sample_and_haplo_class[:, 1],
                     color=color, label=labels.HAPLO_LABELS[haplo_class],
                     alpha=alpha)
    return {'uniq_haplo_classes': uniq_haplo_classes}


def plot_haplo_pop_pcas(fig, pcas_col_width, haplo_pop_pca_locations,
                        aligned_pcoas_df, haplo_classification,
                        pops, legend_pos, common_top_margin, alpha=1):
    num_pca_cols = 3
    num_pca_rows = 3
    pca_col_width = (pcas_col_width / num_pca_cols)
    pca_row_height = 1 / num_pca_rows
    axes_col_widths = [pca_col_width] * num_pca_cols + [1 - pcas_col_width]
    axes_row_heights = [common_top_margin] + [pca_row_height] * num_pca_rows

    for pop, (row_idx, col_idx) in haplo_pop_pca_locations.items():
        axes = matplotlib_support.add_axes(fig, row_idx=row_idx + 1, col_idx=col_idx,
                                            left_margin=0, right_margin=0,
                                            top_margin=0, bottom_margin=0,
                                            axes_col_widths=axes_col_widths,
                                            axes_row_heights=axes_row_heights)
        res = plot_haplo_pca_for_samples(aligned_pcoas_df, haplo_classification,
                                         samples=pops[pop], axes=axes, alpha=alpha)
        axes.set_xlim(HAPLO_PCOAS_X_LIMS)
        axes.set_ylim(HAPLO_PCOAS_Y_LIMS)
        matplotlib_support.set_axes_background(axes)
        matplotlib_support.turn_off_both_axis(axes)
        matplotlib_support.set_axis_color(axes)
        axes.text(HAPLO_PCOAS_X_LIMS[0] + 0.01, sum(HAPLO_PCOAS_Y_LIMS) / 2,
                  labels.LABELS[pop],
                  verticalalignment='center', rotation='vertical',
                  fontsize=20)

    axes = matplotlib_support.add_axes(fig, row_idx=legend_pos[0] + 1, col_idx=legend_pos[1],
                                        left_margin=0, right_margin=0,
                                        top_margin=0, bottom_margin=0,
                                        axes_col_widths=axes_col_widths,
                                        axes_row_heights=axes_row_heights)
    legend_info = [(labels.HAPLO_LABELS[haplo_class], colors.HAPLO_COLORS[haplo_class]) for haplo_class in res['uniq_haplo_classes']]
    legend_info = sorted(legend_info, key=lambda x: x[0])
    labels_, colors_ = zip(*legend_info)
    matplotlib_support.plot_legend(labels_, colors_, axes)
    matplotlib_support.set_axes_background(axes)
    matplotlib_support.turn_off_both_axis(axes)


def plot_hbar_admixtures(admixtures_per_pop, axes,
                         composition_classes_colors):
    y_edges = numpy.arange(admixtures_per_pop.shape[0] + 1)
    y_poss = (y_edges[:-1] + y_edges[1:]) / 2
    height = (y_poss[1] - y_poss[0]) * 0.9

    left = None
    for composition_class in admixtures_per_pop.columns:
        width = admixtures_per_pop.loc[:, composition_class]
        color = composition_classes_colors[composition_class]
        axes.barh(y_poss, width, left=left, align='center', color=color)
        if left is None:
            left = width
        else:
            left += width
    matplotlib_support.set_y_ticks(y_poss, list(admixtures_per_pop.index), axes, fontsize=18)


def plot_structure(axes, prior, k, pops_for_samples, ancestral_pop_colors, pop_order):

    results = parse_faststructure_results(prior)

    admixtures_per_pop = _calc_admixtures_per_pop(results[k]['admixtures'],
                                                  pops_for_samples)
    sorted_index = sorted(admixtures_per_pop.index, key=lambda x: pop_order.index(x))
    admixtures_per_pop = admixtures_per_pop.reindex(sorted_index)
    plot_hbar_admixtures(admixtures_per_pop, axes=axes,
                         composition_classes_colors=ancestral_pop_colors)
    axes.set_xlim((0, 1))
    axes.set_xlabel('Ancestral population freq.', fontsize=18)
    matplotlib_support.set_axes_background(axes)
    matplotlib_support.set_y_ticks_right(axes)

def plot_haplo_compositions(pops, haplo_classification, axes, pop_order):

    haplo_label_mapping = labels.HAPLO_LABELS
    pop_composition_freqs = calc_haplo_pop_composition_freq_dframe(pops, haplo_classification)
    pop_composition_freqs.columns = [haplo_label_mapping[col] for col in pop_composition_freqs.columns]
    sorted_index = sorted(pop_composition_freqs.index, key=lambda x: pop_order.index(labels.LABELS[x]))
    pop_composition_freqs = pop_composition_freqs.reindex(sorted_index)

    haplo_colors = {haplo_label_mapping[haplo]: color for haplo, color in colors.HAPLO_COLORS.items()}

    plot_hbar_admixtures(pop_composition_freqs, axes=axes,
                         composition_classes_colors=haplo_colors)
    axes.set_xlim((0, 1))
    axes.set_xlabel('Haplotype freq.', fontsize=18)
    matplotlib_support.set_axes_background(axes)
    matplotlib_support.turn_off_y_axis(axes)


if __name__ == '__main__':

    out_dir = config.PAPER_FIGURES_DIR
    out_dir.mkdir(exist_ok=True)
    plot_path = out_dir / 'fig3.png'

    pcas_col_width = 0.5
    halpos_pca_alpha = 0.1
    structure_prior = 'simple'
    structure_k = 3

    sample_passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, sample_passports)
    pops_for_samples = {sample: labels.get_long_label(pop) for pop, samples in pops.items() for sample in samples}

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    samples_to_use = {sample for samples in pops.values() for sample in samples}
    res = detected_outliers_and_classify_haplos(variations,
                                                win_params={'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                                                            'win_size': config.HAPLO_WIN_SIZE},
                                                num_wins_to_process=None,
                                                samples_to_use=samples_to_use,
                                                n_dims_to_keep=config.N_DIMS_TO_KEEP,
                                                classification_config=config.CLASSIFICATION_CONFIG,
                                                classification_outlier_config=config.CLASSIFICATION_OUTLIER_CONFIG,
                                                outlier_configs=config.OUTLIER_CONFIGS,
                                                out_dir=config.HAPLO_PCOA_DIR,
                                                pops=pops,
                                                outliers_return_aligned_pcoas=False,
                                                only_outliers=False,
                                                classification_references=config.CLASSIFICATION_REFERENCES,
                                                supervised_classification_config=config.SUPERVISED_CLASSIFICATION_CONFIG,
                                                cache_dir=config.CACHE_DIR)
    aligned_pcoas_df = res['aligned_pcoas_df']
    haplo_classification = res['classification']

    fig = Figure((20, 15))
    FigureCanvas(fig) # Don't remove it or savefig will fail later

    haplo_pop_pca_locations = {'slc_ma': (0, 0), 'sll_mx': (0, 2),
                               'sp_ec': (1, 0), 'slc_ec': (1, 1), 'sll_vint': (1, 2),
                               'sp_pe': (2, 0), 'slc_pe': (2, 1), 'sll_modern': (2, 2)}

    fontsize = 25
    matplotlib_support.write_text_in_figure('Wild', 0.065, 0.95, fig, fontsize=fontsize)
    matplotlib_support.write_text_in_figure('Semi-domesticated', 0.16, 0.95, fig, fontsize=fontsize)
    matplotlib_support.write_text_in_figure('Cultivated', 0.38, 0.95, fig, fontsize=fontsize)

    plot_haplo_pop_pcas(fig, pcas_col_width, haplo_pop_pca_locations,
                        aligned_pcoas_df, haplo_classification,
                        pops, alpha=halpos_pca_alpha, legend_pos=(0, 1),
                        common_top_margin=0.08)

    axes_col_widths = [pcas_col_width, 1 - pcas_col_width]
    axes_row_heights = [0.6, 0.4]
    axes = matplotlib_support.add_axes(fig, row_idx=0, col_idx=1,
                                       left_margin=0, right_margin=0,
                                       top_margin=0, bottom_margin=0,
                                       axes_col_widths=axes_col_widths,
                                       axes_row_heights=axes_row_heights,
                                       projection=ccrs.PlateCarree())
    plot_geo_rank1_for_main_pops(sample_passports, axes=axes)

    axes_col_widths = [0.25, 0.25, 0.22, 0.28]
    structure_pop_order = ['sp_pe', 'sp_x_sp', 'sp_pe_inter-andean',
                           'sp_ec',
                           'slc_co',
                           'slc_ma',
                           'slc_world',
                           'sll_mx', 'sll_vint', 'sll_old_cultivars',
                           'sll_modern',
                           'sp_x_sl',
                           None,
                           'slc_pe',
                           'slc_ec']
    structure_pop_order = list(map(labels.get_long_label, structure_pop_order))
    axes = matplotlib_support.add_axes(fig, row_idx=1, col_idx=3,
                                       left_margin=0, right_margin=0.4,
                                       top_margin=0.01, bottom_margin=0.1,
                                       axes_col_widths=axes_col_widths,
                                       axes_row_heights=axes_row_heights)
    structure_color_mapping = {0: colors.HAPLO_COLORS['sl'],
                               1: colors.HAPLO_COLORS['sp_ecu'],
                               2: colors.HAPLO_COLORS['sp_peru']}
    plot_structure(axes, prior=structure_prior, pops_for_samples=pops_for_samples, k=structure_k,
                   ancestral_pop_colors=structure_color_mapping, pop_order=structure_pop_order)
    matplotlib_support.turn_off_y_axis(axes)
    axes.spines['left'].set_color('#ffffff')

    axes = matplotlib_support.add_axes(fig, row_idx=1, col_idx=2,
                                       left_margin=0.1, right_margin=0.05,
                                       top_margin=0.01, bottom_margin=0.1,
                                       axes_col_widths=axes_col_widths,
                                       axes_row_heights=axes_row_heights)
    plot_haplo_compositions(pops, haplo_classification, axes, pop_order=structure_pop_order)
    axes.spines['left'].set_color('#ffffff')

    fontsize = 25
    matplotlib_support.write_text_in_figure('A', 0.01, 0.97, fig, fontsize=fontsize)
    matplotlib_support.write_text_in_figure('B', 0.55, 0.97, fig, fontsize=fontsize)
    matplotlib_support.write_text_in_figure('C', 0.523, 0.355, fig, fontsize=fontsize)
    matplotlib_support.write_text_in_figure('D', 0.72, 0.355, fig, fontsize=fontsize)

    fig.savefig(plot_path)
