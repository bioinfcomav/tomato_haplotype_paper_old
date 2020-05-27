
import config

from collections import defaultdict

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib
from matplotlib import colors
from matplotlib.cm import ScalarMappable
from matplotlib.offsetbox import AnchoredText

import cartopy.crs as ccrs

from variation.variations import VariationsH5

from plot_geo_map import plot_geo_rank1_for_main_pops
import passport
import plot_thinned_haplos
from pop_building import get_pops
from haplo_pca_plotting import calc_ellipsoids
from colors import PINK_BLUE_CMAP_R2, HAPLO_COLORS, ELLIPSE_COLORS
from haplo_auto_classification import (detected_outliers_and_classify_haplos,
                                       calc_haplo_pop_composition_freq_dframe)
from haplo import get_pop_classification_for_haplos
from haplo_auto_classification_composition_plot import plot_pop_haplo_composition
from faststructure_parse_results import parse_faststructure_results
from faststructure_plot import plot_admixtures, _calc_admixtures_per_pop
from pop_building import get_classifications_for_classification_key_path
import labels


Y_LABEL_TEXT_SIZE = 25
HAPLO_IDS_TEXT_SIZE = 40

HAPLO_PCOAS_X_LIMS = (-0.09, 0.04)
HAPLO_PCOAS_Y_LIMS = (-0.06, 0.04)


def plot_haplo_composition(axes,
                           variations, pops, samples_to_use, pop_order,
                           haplo_label_mapping):
    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}

    out_dir = config.HAPLO_PCOA_DIR
    out_dir.mkdir(exist_ok=True)

    res = detected_outliers_and_classify_haplos(variations,
                                                win_params=win_params,
                                                num_wins_to_process=None,
                                                samples_to_use=samples_to_use,
                                                n_dims_to_keep=config.N_DIMS_TO_KEEP,
                                                classification_config=config.CLASSIFICATION_CONFIG,
                                                classification_outlier_config=config.CLASSIFICATION_OUTLIER_CONFIG,
                                                outlier_configs=config.OUTLIER_CONFIGS,
                                                out_dir=out_dir,
                                                classification_references=config.CLASSIFICATION_REFERENCES,
                                                supervised_classification_config=config.SUPERVISED_CLASSIFICATION_CONFIG,
                                                cache_dir=cache_dir)
    haplo_classification = res['classification']
    haplo_pop_classification = get_pop_classification_for_haplos(res['aligned_pcoas_df'].index, pops)
    pop_composition_freqs = calc_haplo_pop_composition_freq_dframe(pops, haplo_classification)
    pop_composition_freqs.columns = [haplo_label_mapping[col] for col in pop_composition_freqs.columns]

    haplo_colors = {haplo_label_mapping[haplo]: color for haplo, color in HAPLO_COLORS.items()}

    plot_admixtures(pop_composition_freqs, axes=axes, composition_classes_colors=haplo_colors,
                    sample_order=pop_order)


def plot_structure(axes, prior, k, pops_for_samples, ancestral_pop_colors, pop_order):

    results = parse_faststructure_results(prior)

    admixtures_per_pop = _calc_admixtures_per_pop(results[k]['admixtures'],
                                                  pops_for_samples)
    plot_admixtures(admixtures_per_pop, axes=axes, composition_classes_colors=ancestral_pop_colors,
                    sample_order=pop_order)


def plot_figure1(plot_path, sample_passports,
                 variations, dist_threshold, samples_to_use, cache_dir,
                 structure_prior, structure_k,
                 pops_for_haplo_classification,
                 haplo_label_mapping,
                 ):
    fig = Figure((20, 20))
    FigureCanvas(fig) # Don't remove it or savefig will fail later

    axes_id_letters_text_size = 40
    nrows = 3
    ncols = 2
    height_ratios = [0.9, 0.1, 1]

    grid_spec = matplotlib.gridspec.GridSpec(nrows, ncols, figure=fig,
                                             height_ratios=height_ratios)

    map_axes = fig.add_subplot(grid_spec[2, 0], projection=ccrs.PlateCarree(), zorder=1)
    plot_geo_rank1_for_main_pops(sample_passports, axes=map_axes)

    text = AnchoredText('C', prop=dict(size=axes_id_letters_text_size), loc='upper left', frameon=False)
    map_axes.add_artist(text)
    text.set_zorder(100)

    if False:
        haplos_hist_axes = fig.add_subplot(grid_spec[0, 0])
        res = plot_thinned_haplos.plot_thinned_haplos(variations, axes=haplos_hist_axes,
                                                    dist_threshold=dist_threshold,
                                                    samples_to_use=samples_to_use,
                                                    pops=pops_for_haplo_classification,
                                                    cache_dir=cache_dir)
        haplos_hist_axes.set_facecolor('white')
        haplos_hist_axes.set_xlabel('Dim. 1')
        haplos_hist_axes.set_ylabel('Dim. 2')
        haplos_hist_axes.text(0.0075, 0.0070, 'SL', fontsize=HAPLO_IDS_TEXT_SIZE, color='#782d50', zorder=90)
        haplos_hist_axes.text(-0.050, -0.045, 'SP EC', fontsize=HAPLO_IDS_TEXT_SIZE, color='#1f6e8c', zorder=90)
        haplos_hist_axes.text(-0.065, -0.02, 'SP PE', fontsize=HAPLO_IDS_TEXT_SIZE, color='#327795', zorder=90)
        colorbar_axes = fig.add_subplot(grid_spec[1, 0])
        fig.colorbar(res['hist2d_result'][3], cax=colorbar_axes, orientation='horizontal')
    else:
        haplos_hist_axes = fig.add_subplot(grid_spec[0:2, 0])
        res = plot_thinned_haplos.plot_thinned_haplos2(variations, axes=haplos_hist_axes,
                                                    dist_threshold=dist_threshold,
                                                    samples_to_use=samples_to_use,
                                                    pops=pops_for_haplo_classification,
                                                    cache_dir=cache_dir,
                                                    alpha=0.05,
                                                    haplo_colors=HAPLO_COLORS)
        haplos_hist_axes.set_xlim(HAPLO_PCOAS_X_LIMS)
        haplos_hist_axes.set_ylim(HAPLO_PCOAS_Y_LIMS)
        haplos_hist_axes.set_facecolor('white')

    text = AnchoredText('A', prop=dict(size=axes_id_letters_text_size), loc='upper left', frameon=False)
    haplos_hist_axes.add_artist(text)
    text.set_zorder(100)

    pops_for_samples = get_classifications_for_classification_key_path(sample_passports, config.RANK1)
    pops_for_samples = {sample: labels.get_long_label(pop) for sample, pop in pops_for_samples.items()}
    pops = defaultdict(list)
    for sample_name, pop in pops_for_samples.items():
        pops[pop].append(sample_name)

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
                
    structure_color_mapping = {0: HAPLO_COLORS['sl'], 1: HAPLO_COLORS['sp_ecu'], 2: HAPLO_COLORS['sp_peru']}

    structure_axes = fig.add_subplot(grid_spec[2, 1])
    plot_structure(structure_axes, prior=structure_prior, pops_for_samples=pops_for_samples, k=structure_k,
                   ancestral_pop_colors=structure_color_mapping, pop_order=structure_pop_order)
    structure_axes.tick_params(axis='x', labelsize=20)
    structure_axes.set_ylabel('Ancestral population freq.', fontsize=Y_LABEL_TEXT_SIZE)
    structure_axes.yaxis.set_label_position('right')
    structure_axes.tick_params(axis='y',
                                       which='both',
                                       right=True,
                                       left=False, 
                                       labelleft=False,
                                       labelright=True
                                       )
    text = AnchoredText('D', prop=dict(size=axes_id_letters_text_size), loc='upper left', frameon=False)
    structure_axes.add_artist(text)
    text.set_zorder(100)

    haplo_composition_axes = fig.add_subplot(grid_spec[0:2, 1])
    plot_haplo_composition(haplo_composition_axes,
                           variations, pops, samples_to_use, pop_order=structure_pop_order,
                           haplo_label_mapping=haplo_label_mapping)
    haplo_composition_axes.tick_params(axis='x',          # changes apply to the x-axis
                                       which='both',      # both major and minor ticks are affected
                                       bottom=False,      # ticks along the bottom edge are off
                                       top=False,         # ticks along the top edge are off
                                       labelbottom=False)
    haplo_composition_axes.set_ylabel('Haplotype freq.', fontsize=Y_LABEL_TEXT_SIZE)
    haplo_composition_axes.yaxis.set_label_position('right')
    haplo_composition_axes.tick_params(axis='y',
                                       which='both',
                                       right=True,
                                       left=False, 
                                       labelleft=False,
                                       labelright=True
                                       )
    haplo_composition_axes.legend(prop={'size': 20})
    text = AnchoredText('B', prop=dict(size=axes_id_letters_text_size), loc='upper left', frameon=False)
    haplo_composition_axes.add_artist(text)
    text.set_zorder(100)

    fig.savefig(plot_path)


if __name__ == '__main__':

    cache_dir = config.CACHE_DIR
    dist_threshold = config.CLASSIFICATION_CONFIG['thinning_dist_threshold']

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    imputed_variations = VariationsH5(str(vars_path), 'r')

    out_dir = config.PAPER_FIGURES_DIR
    out_dir.mkdir(exist_ok=True)
    plot_path = out_dir / 'fig1.svg'

    sample_passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, sample_passports)
    samples_to_use = {sample for samples in pops.values() for sample in samples}

    plot_figure1(plot_path, sample_passports=sample_passports,
                 variations=imputed_variations, dist_threshold=dist_threshold,
                 samples_to_use=samples_to_use, cache_dir=cache_dir,
                 pops_for_haplo_classification=pops,
                 structure_prior='simple', structure_k=3,
                 haplo_label_mapping=labels.HAPLO_LABELS)
