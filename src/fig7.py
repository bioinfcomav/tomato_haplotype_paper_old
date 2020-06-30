
import config

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5

import passport
import pop_building
from haplo_auto_classification import detected_outliers_and_classify_haplos
import colors
import labels
from haplo_priv import (count_num_shared_and_priv_haplos,
                        plot_shared_and_priv_haplo_counts)
import matplotlib_support


if __name__ == '__main__':

    debug = False
    classify_haplos = True
    if debug:
        num_wins_to_process = 2
    else:
        num_wins_to_process = None

    win_params={'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                'win_size': config.HAPLO_WIN_SIZE}

    sample_passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, sample_passports)

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    if classify_haplos:
        res = detected_outliers_and_classify_haplos(variations,
                                                    win_params=win_params,
                                                    num_wins_to_process=None,
                                                    samples_to_use=variations.samples,
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
        haplo_classification = res['classification']
    else:
        haplo_classification = None

    pop1 = 'slc_ma'
    pop2 = 'slc_pe'
    pop3 = 'sll_mx'
    pop_names = pop1, pop2, pop3
    pop_combination_order = [('slc_ma', 'slc_pe'), ('slc_pe', 'sll_mx'), ('slc_ma', 'sll_mx'),
                              'slc_ma', 'slc_pe', 'sll_mx']
    color_schema = colors.ColorSchema(colors.HAPLO_COLORS)

    pops = {pop: pops[pop] for pop in pop_names}

    counts = count_num_shared_and_priv_haplos(pops,
                                              variations,
                                              win_params=win_params,
                                              haplo_classification=haplo_classification,
                                              num_wins_to_process=num_wins_to_process)

    out_dir = config.FIGURES_DIR

    plot_path = out_dir / f'fig7.svg'
    fig = Figure((10, 7))
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)
    res = plot_shared_and_priv_haplo_counts(counts, axes, sorted_pop_names=pop_names, ignore_shared_by_all_pops=True, pop_combination_order=pop_combination_order)
    matplotlib_support.plot_legend([labels.HAPLO_LABELS[klass] for klass in res['haplo_classes']],
                                   [color_schema[klass] for klass in res['haplo_classes']], axes)
    matplotlib_support.set_axes_background(axes)
    fig.tight_layout()
    fig.savefig(plot_path)
