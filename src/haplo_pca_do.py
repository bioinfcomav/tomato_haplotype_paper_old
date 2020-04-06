
import config

from collections import defaultdict

import numpy
import pandas

from variation.variations import VariationsH5

from passport import get_sample_passports
from haplo_pca import do_pcoas_along_the_genome
from procrustes import align_pcas_using_procrustes
from pop_building import get_pops
from haplo_pca_plotting import (plot_haplo_pcas, write_pcas_curly_file,
                                plot_pcas_per_pop, plot_pcas_per_sample)


def get_sample_selection_criteria():

    rank1 = config.RANK1

    criteria = []
    samples_to_remove = []
    samples_to_keep = []

    return {'criteria': criteria, 'samples_to_remove': samples_to_remove,
            'samples_to_keep': samples_to_keep}


def get_haplotypes_to_exclude(path):
    haplos_to_exclude = defaultdict(list)

    if not path.exists():
        return {}

    for line in path.open('rt'):
        if line.startswith('Label'):
            continue
        line = line.strip()
        if not line:
            continue
        items = line.strip().split()[0].split('%')
        chrom = f'SL2.50ch{items[0]}'
        win_start = int(items[1])
        key = chrom, win_start
        sample = items[2]
        haploid_idx = int(items[3])
        haplos_to_exclude[key].append((sample, haploid_idx))
    return haplos_to_exclude


def get_haplotypes_to_include(path):
    return get_haplotypes_to_exclude(path)


if __name__ == '__main__':

    debug = True

    if debug:
        num_wins_to_process = 2
    else:
        num_wins_to_process = None

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

    pcoas = do_pcoas_along_the_genome(variations, win_params,
                                      num_wins_to_process=num_wins_to_process,
                                      samples=samples_to_use, n_dims_to_keep=3)

    aligned_pcoas = list(align_pcas_using_procrustes(pcoas))

    out_dir = config.HAPLO_PCOA_DIR
    out_dir.mkdir(exist_ok=True)

    haplotypes_to_include_path = out_dir / 'haplotypes_to_include.txt'
    haplotypes_to_include = get_haplotypes_to_include(haplotypes_to_include_path)

    haplotypes_to_exclude_path = out_dir / 'haplotypes_to_exclude.txt'
    haplotypes_to_exclude = get_haplotypes_to_exclude(haplotypes_to_exclude_path)

    haplo_classification = None

    write_pcas_curly_file(aligned_pcoas, out_dir, pops, haplo_classification)

    ellipsoids = None

    res = plot_haplo_pcas(aligned_pcoas, out_dir, populations=pops,
                          ellipsoids=ellipsoids)

    per_pop_out_dir = out_dir / 'per_pop'
    per_pop_out_dir.mkdir(exist_ok=True)
    plot_pcas_per_pop(aligned_pcoas, per_pop_out_dir, res['x_lims'], res['y_lims'],
                      color_schema=res['color_schema'],
                      populations=pops, ellipsoids=ellipsoids)

    per_sample_out_dir = out_dir / 'per_sample'
    per_sample_out_dir.mkdir(exist_ok=True)
    plot_pcas_per_sample(aligned_pcoas, per_sample_out_dir,
                         res['x_lims'], res['y_lims'],
                         populations=pops,
                         ellipsoids=ellipsoids)
