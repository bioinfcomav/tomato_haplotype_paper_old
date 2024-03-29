
import config

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5

from passport import get_sample_passports
from pop_building import get_pops
from haplo_auto_classification import (detected_outliers_and_classify_haplos,
                                       calc_haplo_pop_composition_freq)
from haplo import get_pop_classification_for_haplos, parse_haplo_id


def plot_pop_haplo_composition(pops, haplo_classification, haplo_pop_classification,
                               plot_path=None, axes=None):

    pop_composition_freqs = calc_haplo_pop_composition_freq(pops, haplo_classification)

    haplo_classes = sorted({klass for compositions in pop_composition_freqs.values() for klass in compositions},
                           key=str)

    pop_names = sorted(pops.keys(), key=str)

    if axes is None:
        fig = Figure()
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)
        savefig = True
    else:
        savefig = False

    x_poss = numpy.arange(len(pop_names))
    width = x_poss[1] - x_poss[0]
    bottoms = None
    for klass in haplo_classes:
        pop_freqs = numpy.array([pop_composition_freqs[pop][klass] for pop in pop_names])
        axes.bar(x_poss, pop_freqs, width, bottoms, label=klass)
        if bottoms is None:
            bottoms = pop_freqs
        else:
            bottoms += pop_freqs

    axes.set_xticklabels(pop_names, rotation=45, ha='right')
    axes.set_xticks(x_poss)

    axes.legend()

    if savefig:
        fig.tight_layout()
        fig.savefig(str(plot_path))


if __name__ == '__main__':

    debug = False

    num_wins_to_process = None
    cache_dir = config.CACHE_DIR
    only_outliers = False
    outliers_return_aligned_pcoas = False

    win_params = {'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                  'win_size': config.HAPLO_WIN_SIZE}

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = get_sample_passports()

    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, passports)

    samples_to_use = {sample for samples in pops.values() for sample in samples}
    if config.SKIP_SAMPLES_WITH_NO_GENOTYPE:
        samples_to_use = samples_to_use.intersection(variations.samples)
    samples_to_use = sorted(samples_to_use)

    out_dir = config.HAPLO_PCOA_DIR
    out_dir.mkdir(exist_ok=True)

    res = detected_outliers_and_classify_haplos(variations,
                                                win_params=win_params,
                                                num_wins_to_process=num_wins_to_process,
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

    plot_path = out_dir / 'pop_haplo_composition.svg'
    plot_pop_haplo_composition(pops, haplo_classification, haplo_pop_classification, plot_path)
