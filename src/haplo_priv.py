
import config

from collections import defaultdict, Counter

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5

from haplo import generate_uniq_haplos_along_genome, create_haplo_id
import passport
import pop_building
from haplo_auto_classification import detected_outliers_and_classify_haplos
import colors
import matplotlib_support
import fig_style
import labels


def count_num_shared_and_priv_haplos(pops, variations, win_params, haplo_classification=None,
                                     num_wins_to_process=None):

    pop_names = sorted(pops.keys())
    samples_in_pops = {pop: set(samples) for pop, samples in pops.items()}

    if haplo_classification is None:
        haplo_classification = {}

    counts = defaultdict(Counter)
    for haplo_info in generate_uniq_haplos_along_genome(variations, win_params, num_wins_to_process=num_wins_to_process):
        win_start = haplo_info['win_start']
        chrom = haplo_info['chrom']
        uniq_haplo_ids_that_correpond_to_all_sample_haplos = haplo_info['uniq_haplo_ids_that_correpond_to_all_haplos']

        haplos_in_each_pop = {pop: set() for pop in pop_names}
        haplo_classes = {}
        for haplo_id, (sample, haploid_idx) in zip(uniq_haplo_ids_that_correpond_to_all_sample_haplos.values,
                                                   uniq_haplo_ids_that_correpond_to_all_sample_haplos.index):
            for pop, samples_in_pop in samples_in_pops.items():
                if sample in samples_in_pop:
                    haplos_in_each_pop[pop].add(haplo_id)
            haplo_classes[haplo_id] = haplo_classification.get(create_haplo_id(chrom, win_start, sample, haploid_idx), None)

        for haplo_id in uniq_haplo_ids_that_correpond_to_all_sample_haplos.values:
            pops_in_which_haplo_is_found = tuple(pop for pop in pop_names if haplo_id in haplos_in_each_pop[pop])
            if not pops_in_which_haplo_is_found:
                continue
            haplo_class = haplo_classes[haplo_id]
            counts[pops_in_which_haplo_is_found][haplo_class] += 1
    return counts


def plot_shared_and_priv_haplo_counts(counts, axes, sorted_pop_names=None, haplo_color_schema=None,
                                      ignore_shared_by_all_pops=False, pop_combination_order=None):

    pops = {pop for pops in counts.keys() for pop in pops}

    if sorted_pop_names is None:
        sorted_pop_names = sorted(pops)

    assert len(pops.intersection(sorted_pop_names)) == len(sorted_pop_names)

    haplo_classes = sorted({haplo_class for counts_in_pops in counts.values() for haplo_class in counts_in_pops.keys()})

    if haplo_color_schema is None:
        haplo_color_schema = colors.ColorSchema(colors.HAPLO_COLORS)

    x_pos_edges = numpy.arange(len(sorted_pop_names) + 1)
    pops_x_pos = (x_pos_edges[:-1] + x_pos_edges[1:]) / 2
    left_pops_x_poss = dict(zip(sorted_pop_names, x_pos_edges[:-1]))

    width = x_pos_edges[1] - x_pos_edges[0]

    pop_combinations = list(reversed(sorted(counts.keys(), key=lambda x: len(x))))

    if ignore_shared_by_all_pops:
        pop_combinations = [pop_combination for pop_combination in pop_combinations if len(pop_combination) < len(sorted_pop_names)]

    if pop_combination_order:
        sort_order = {pop_combination: idx for idx, pop_combination in enumerate(reversed(pop_combination_order))}
        pop_combinations = list(reversed(sorted(pop_combinations, key=lambda x: sort_order.get(x, -1))))

    accumulated_haplo_counts = 0
    for pop_combination in pop_combinations:
        counts_for_each_haplo_class = counts[pop_combination]
        for haplo_class in haplo_classes:
            try:
                haplo_counts = counts_for_each_haplo_class[haplo_class]
            except KeyError:
                continue
            bottom = accumulated_haplo_counts
            height = haplo_counts
            color = haplo_color_schema[haplo_class]

            for pop in pop_combination:
                left_x_pos = left_pops_x_poss[pop]
                axes.bar(left_x_pos, height, width=width, bottom=bottom, align='edge', color=color)
            accumulated_haplo_counts += haplo_counts

    axes.set_xlim((x_pos_edges[0], x_pos_edges[-1]))
    axes.set_ylim((0, accumulated_haplo_counts))
    axes.set_ylabel('Haplotypes', fontsize=fig_style.X_LABEL_SMALL_FONT_SIZE)

    x_labels = [labels.LABELS[pop] for pop in sorted_pop_names]
    matplotlib_support.set_x_ticks(pops_x_pos, x_labels, axes, rotation=45, fontsize=fig_style.X_LABEL_SMALL_FONT_SIZE)
    return {'haplo_classes': haplo_classes}


if __name__ == '__main__':

    debug = True
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

    out_dir = config.TMP_DIR
    pops_str = '-'.join(pop_names)
    plot_path = out_dir / f'haplo_counts_for_{pops_str}.svg'
    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)
    plot_shared_and_priv_haplo_counts(counts, axes, sorted_pop_names=pop_names, pop_combination_order=pop_combination_order)
    matplotlib_support.set_axes_background(axes)
    fig.tight_layout()
    fig.savefig(plot_path)

    plot_path = out_dir / f'haplo_counts_for_{pops_str}.ignore_common.svg'
    fig = Figure((10, 7))
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)
    res = plot_shared_and_priv_haplo_counts(counts, axes, sorted_pop_names=pop_names, ignore_shared_by_all_pops=True, pop_combination_order=pop_combination_order)
    matplotlib_support.plot_legend([labels.HAPLO_LABELS[klass] for klass in res['haplo_classes']],
                                   [color_schema[klass] for klass in res['haplo_classes']], axes)
    matplotlib_support.set_axes_background(axes)
    fig.tight_layout()
    fig.savefig(plot_path)
