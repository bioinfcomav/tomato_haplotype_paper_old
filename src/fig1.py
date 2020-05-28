

import config

import numpy

from variation.variations import VariationsH5

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import passport
from pop_building import get_pops
import plot_thinned_haplos
from colors import HAPLO_COLORS
import labels
import matplotlib_support
from haplo_pca_plotting import filter_aligned_pcoas_df_for_samples
from haplo_auto_classification import detected_outliers_and_classify_haplos

HAPLO_PCOAS_X_LIMS = (-0.09, 0.04)
HAPLO_PCOAS_Y_LIMS = (-0.06, 0.04)
LEGEND_FONT_SIZE = 20
X_LABEL_SMALL_FONT_SIZE = 20
LEGEND_MARKER_SIZE = 15


def plot_haplo_pca(axes, sample_passports,
                   variations, dist_threshold, samples_to_use, cache_dir,
                   structure_prior, structure_k,
                   pops_for_haplo_classification,
                   haplo_label_mapping):

    res = plot_thinned_haplos.plot_thinned_haplos2(variations, axes=haplos_pca_axes,
                                                   dist_threshold=dist_threshold,
                                                   samples_to_use=samples_to_use,
                                                   pops=pops_for_haplo_classification,
                                                   cache_dir=cache_dir,
                                                   alpha=0.05,
                                                   haplo_colors=HAPLO_COLORS)
    haplos_pca_axes.set_xlim(HAPLO_PCOAS_X_LIMS)
    haplos_pca_axes.set_ylim(HAPLO_PCOAS_Y_LIMS)
    haplos_pca_axes.set_facecolor('white')

    legend_elements = []
    for haplo_class in res['haplo_classes']:
        color = res['colors'][haplo_class]
        element = Line2D([0], [0], marker='o', color=color,
                         label=labels.HAPLO_LABELS[haplo_class],
                         markersize=LEGEND_MARKER_SIZE, linestyle='')
        legend_elements.append(element)
    legend = haplos_pca_axes.legend(handles=legend_elements, prop={'size': LEGEND_FONT_SIZE})
    frame = legend.get_frame()
    frame.set_color('#ffffff')
    frame.set_edgecolor('#cccccc')
    frame.set_alpha(0.8)

    haplos_pca_axes.set_xlabel('Dim. 1', fontsize=X_LABEL_SMALL_FONT_SIZE)
    haplos_pca_axes.set_ylabel('Dim. 2', fontsize=X_LABEL_SMALL_FONT_SIZE)

    fig.savefig(plot_path)


def plot_haplo_pca_for_sample(axes, variations, sample, samples_to_use, pops, cache_dir):

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
                                                cache_dir=cache_dir)

    aligned_pcoas_df_for_sample = filter_aligned_pcoas_df_for_samples(res['aligned_pcoas_df'],
                                                                      [sample])
    
    haplo_classification = res['classification']
    haplo_classes = numpy.array([haplo_classification[haplo_id] for haplo_id in aligned_pcoas_df_for_sample.index])
    uniq_haplo_classes = sorted(set(haplo_classes))
    
    for haplo_class in uniq_haplo_classes:
        mask = haplo_classes == haplo_class
        aligned_pcoas_for_sample_and_haplo_class = aligned_pcoas_df_for_sample.values[mask, :]
        color = HAPLO_COLORS.get(haplo_class, None)
        axes.scatter(aligned_pcoas_for_sample_and_haplo_class[:, 0],
                     aligned_pcoas_for_sample_and_haplo_class[:, 1],
                     color=color, label=labels.HAPLO_LABELS[haplo_class])
    axes.set_xlabel(f'Dim 1', fontsize=X_LABEL_SMALL_FONT_SIZE)
    axes.set_ylabel(f'Dim 2', fontsize=X_LABEL_SMALL_FONT_SIZE)
    axes.set_facecolor('white')


if __name__ == '__main__':

    cache_dir = config.CACHE_DIR
    dist_threshold = config.CLASSIFICATION_CONFIG['thinning_dist_threshold']
    sample_for_haplo_pca = 'cervil'

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    imputed_variations = VariationsH5(str(vars_path), 'r')

    out_dir = config.PAPER_FIGURES_DIR
    out_dir.mkdir(exist_ok=True)
    plot_path = out_dir / 'fig1.svg'

    sample_passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, sample_passports)
    samples_to_use = {sample for samples in pops.values() for sample in samples}

    fig = Figure((10, 20))
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes_row_heights = [0.5, 0.5]

    haplos_pca_axes = matplotlib_support.add_axes(fig, row_idx=0, axes_row_heights=axes_row_heights)

    plot_haplo_pca(haplos_pca_axes, sample_passports=sample_passports,
                   variations=imputed_variations, dist_threshold=dist_threshold,
                   samples_to_use=samples_to_use, cache_dir=cache_dir,
                   pops_for_haplo_classification=pops,
                   structure_prior='simple', structure_k=3,
                   haplo_label_mapping=labels.HAPLO_LABELS)

    axes2 = matplotlib_support.add_axes(fig, row_idx=1, axes_row_heights=axes_row_heights)
    plot_haplo_pca_for_sample(axes2, imputed_variations, sample_for_haplo_pca, samples_to_use, pops, cache_dir)

    fig.savefig(plot_path)
