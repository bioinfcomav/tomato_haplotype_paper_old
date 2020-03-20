
import config

from collections import defaultdict

import numpy
import pandas

from variation.variations import VariationsH5

from passport import get_sample_passports
from haplo_pca import do_pcoas_along_the_genome
from procrustes import align_pcas_using_procrustes
from pca import write_curlywhirly_file
from pop_building import get_pops


def write_pcas_curly_file(pcoas, out_dir, populations=None, haplo_classification=None):

    pop_for_samples = {sample: pop for pop, samples in populations.items() for sample in samples}

    path = out_dir / 'pcoas_along_the_genome.curly'

    sample_names = []
    pops = []
    all_projections = []
    for pcoa_idx, pcoa in enumerate(pcoas):
        projections = pcoa['projections']
        chrom = pcoa['chrom'][-2:]

        if populations:
            sample_names.extend([f'{chrom}%{pcoa["win_start"]}%{sample}%{haploid_idx}' for sample, haploid_idx in projections.index])
            pops.extend([pop_for_samples[sample] for sample, _ in projections.index])

        all_projections.append(projections.values)

    all_projections = numpy.vstack(all_projections)

    categories = {}

    if populations:
        pop_classification = dict(zip(sample_names, pops))
        categories['pop'] = pop_classification

    if haplo_classification:
        haplo_classification = {haplo_id: pop for pop, haplo_ids in haplo_classification.items() for haplo_id in haplo_ids}
        categories['haplo_class'] = haplo_classification

    if not categories:
        categories = None

    all_projections = pandas.DataFrame(all_projections, index=sample_names)
    write_curlywhirly_file(all_projections, path,
                           categories=categories)


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

    vars_path = config.TIER1_PHASED_AND_IMPUTED_LOW_QUAL_09_MISSING_085
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

    aligned_pcoas = align_pcas_using_procrustes(pcoas)

    out_dir = config.HAPLO_PCOA_DIR
    out_dir.mkdir(exist_ok=True)

    haplotypes_to_include_path = out_dir / 'haplotypes_to_include.txt'
    haplotypes_to_include = get_haplotypes_to_include(haplotypes_to_include_path)

    haplotypes_to_exclude_path = out_dir / 'haplotypes_to_exclude.txt'
    haplotypes_to_exclude = get_haplotypes_to_exclude(haplotypes_to_exclude_path)

    haplo_classification = None

    write_pcas_curly_file(aligned_pcoas, out_dir, pops, haplo_classification)

    # plot_pcas(aligned_pcoas, out_dir, populations, ellipsoids)
