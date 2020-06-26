
import config

from collections import OrderedDict

import numpy

from variation.variations import VariationsH5

import passport
import pop_building
import genome
import introgressions
from introgressions_and_kegg import _calc_gene_introgression_freqs


def calc_introgression_freqs_for_genes(variations, introgression_freqs, genes, method='mean_highest'):

    vars_index = variations.pos_index

    mean_gene_introgression_freqs = {}
    for gene_id in genes.gene_ids:
        try:
            gene_introgression_freqs = _calc_gene_introgression_freqs(gene_id, genes, vars_index, introgression_freqs)
        except KeyError:
            continue

        if method == 'mean_highest':
            q90 = numpy.nanpercentile(gene_introgression_freqs, 90)
            vars_with_most_introgression = gene_introgression_freqs[gene_introgression_freqs > q90]
            if vars_with_most_introgression.size:
                mean_gene_introgression_freqs[gene_id] = numpy.nanmean(vars_with_most_introgression)

    return mean_gene_introgression_freqs 


def get_genes_with_most_introgressions(variations, introgression_freqs, genes, num_genes=10):
    introgression_freqs = calc_introgression_freqs_for_genes(variations, introgression_freqs, genes)
    gene_ids = list(reversed(sorted(introgression_freqs.keys(), key=lambda gene_id: introgression_freqs[gene_id])))[:num_genes]

    gene_introgression_freqs = OrderedDict([(gene_id, introgression_freqs[gene_id]) for gene_id in gene_ids])
    return gene_introgression_freqs


if __name__ == '__main__':
    vars_path = config.WORKING_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops_rank1 = pop_building.get_pops(pops_descriptions, passports)

    pops_descriptions = {config.RANK2: config.ALL_POPS}
    pops_rank2 = pop_building.get_pops(pops_descriptions, passports)

    win_size = 10000
    method = 'mean_highest'

    genes = genome.Genes()

    founder_pop = 'slc_ma'

    target_and_introgression_source = [('slc_pe', config.RANK1, 'sp_pe'),
                                       ('slc_ec', config.RANK1, 'sp_ec'),
                                       ]
    out_dir = config.FIGURES_DIR
    out_dir.mkdir(exist_ok=True)

    for target_pop, target_pop_rank, introgression_source_pop in target_and_introgression_source:
        print(f'Introgressions in {target_pop}')
        if target_pop_rank == config.RANK1:
            samples_in_pop_with_possible_introgressions = pops_rank1[target_pop]
        else:
            samples_in_pop_with_possible_introgressions = pops_rank2[target_pop]

        samples_in_founder_pop = pops_rank1[founder_pop]
        samples_in_introgression_source_pop = pops_rank1[introgression_source_pop]

        introgession_config = {'samples_in_pop_with_introgressions': samples_in_pop_with_possible_introgressions,
                               'samples_in_founder_pop': samples_in_founder_pop,
                               'samples_in_introgression_source_pop': samples_in_introgression_source_pop,
                               'freq_threshold_to_consider_allele_present_in_founder_pop' : config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_PRESENT,
                               'freq_threshold_to_consider_allele_common_in_introgression_source_pop': config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_COMMON,
                               }
        introgression_freqs = introgressions.calc_introgession_freq_for_vars(variations, introgession_config)

        gene_introgression_freqs = get_genes_with_most_introgressions(variations, introgression_freqs, genes, num_genes=5)
        print(gene_introgression_freqs)
