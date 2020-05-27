
import config

from variation.variations import VariationsH5

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import passport
from pop_building import get_pops
import plot_thinned_haplos
from colors import HAPLO_COLORS
import labels


HAPLO_PCOAS_X_LIMS = (-0.09, 0.04)
HAPLO_PCOAS_Y_LIMS = (-0.06, 0.04)
LEGEND_FONT_SIZE = 20
X_LABEL_SMALL_FONT_SIZE = 20


def plot_figure1(plot_path, sample_passports,
                 variations, dist_threshold, samples_to_use, cache_dir,
                 structure_prior, structure_k,
                 pops_for_haplo_classification,
                 haplo_label_mapping,
                 ):
    fig = Figure((10, 10))
    FigureCanvas(fig) # Don't remove it or savefig will fail later

    haplos_pca_axes = fig.add_subplot(111)

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
                         markersize=15, linestyle='')
        legend_elements.append(element)
    legend = haplos_pca_axes.legend(handles=legend_elements, prop={'size': LEGEND_FONT_SIZE})
    frame = legend.get_frame()
    frame.set_color('#ffffff')
    frame.set_edgecolor('#cccccc')
    frame.set_alpha(0.8)

    haplos_pca_axes.set_xlabel('Dim. 1', fontsize=X_LABEL_SMALL_FONT_SIZE)
    haplos_pca_axes.set_ylabel('Dim. 2', fontsize=X_LABEL_SMALL_FONT_SIZE)

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
