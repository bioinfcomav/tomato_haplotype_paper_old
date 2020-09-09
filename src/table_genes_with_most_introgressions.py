
import config

import csv
from collections import OrderedDict

from variation.variations import VariationsH5

import passport
import pop_building
import genome
import introgressions
from introgressions_for_genes import get_genes_with_most_introgressions
import kegg
import labels


def get_genes_with_most_introgressions_for_target_pops(variations, target_and_introgression_source, pops,
                                                       upstream_region=0, cache_dir=None):

    genes_with_most_introgessions = {}
    for target_pop, introgression_source_pop in target_and_introgression_source:
        print(f'Introgressions in {target_pop}')
        samples_in_pop_with_possible_introgressions = pops[target_pop]
        samples_in_founder_pop = pops[founder_pop]
        samples_in_introgression_source_pop = pops[introgression_source_pop]

        introgession_config = {'samples_in_pop_with_introgressions': sorted(samples_in_pop_with_possible_introgressions),
                               'samples_in_founder_pop': sorted(samples_in_founder_pop),
                               'samples_in_introgression_source_pop': sorted(samples_in_introgression_source_pop),
                               'freq_threshold_to_consider_allele_present_in_founder_pop' : config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_PRESENT,
                               'freq_threshold_to_consider_allele_common_in_introgression_source_pop': config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_COMMON,
                               }

        gene_introgression_freqs = get_genes_with_most_introgressions(variations,
                                                                      introgession_config=introgession_config,
                                                                      genes=genes, num_genes=100,
                                                                      upstream_region=upstream_region,
                                                                      cache_dir=cache_dir)
        genes_with_most_introgessions[target_pop] = gene_introgression_freqs
    return genes_with_most_introgessions


def remove_suffix(string, suffix):
    if string.endswith(suffix):
        return string[:-len(suffix)]

    raise ValueError(f'suffix {suffix} not in {string}')


def write_genes_with_introgressions_table(fhand, genes_with_most_introgessions, genes):

    writer = csv.writer(fhand, delimiter='\t')
    writer.writerow(['Population', 'Gene id', 'Chromosome', 'Start', 'End', 'Introgression freq.', 'Annotated function', 'Kegg pathways'])

    kegg_genes = kegg.KeggGenes(cache_dir=config.CACHE_DIR)
    for pop in sorted(genes_with_most_introgessions.keys()):
        for gene_id, introgression_freq in genes_with_most_introgessions[pop].items():
            gene = genes.get_gene(gene_id)
            function = genes.get_annotated_function(gene_id)

            try:
                kegg_pathways = kegg_genes.get_patways_for_solcap_id(gene_id)
                kegg_pathways = [remove_suffix(pathway['description'], ' - Solanum lycopersicum (tomato)') for pathway in kegg_pathways]
                kegg_pathways = ','.join(kegg_pathways)
            except KeyError:
                kegg_pathways = ''
            #print(pop, gene_id, introgression_freq, function, kegg_pathways)
            writer.writerow([labels.LABELS[pop], gene_id,
                             gene['Chromosome'], str(gene['Start']), str(gene['End']),
                             introgression_freq, function, kegg_pathways])

if __name__ == '__main__':
    vars_path = config.TIER1_H5_LOWQ_085

    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, passports)

    win_size = 10000
    method = 'mean_highest'

    genes = genome.Genes()

    founder_pop = 'slc_ma'

    target_and_introgression_source = [('slc_pe', 'sp_pe'),
                                       ('slc_ec', 'sp_ec'),
                                       ]
    out_dir = config.FIGURES_DIR
    out_dir.mkdir(exist_ok=True)

    genes_with_most_introgessions = get_genes_with_most_introgressions_for_target_pops(variations,
                                                                                       target_and_introgression_source,
                                                                                       pops,
                                                                                       upstream_region=config.UPSTREAM_REGION_IN_BP,
                                                                                       cache_dir=config.CACHE_DIR)

    path = config.TABLE_GENES_WITH_MANY_INTROGRESSIONS
    write_genes_with_introgressions_table(path.open('wt'), genes_with_most_introgessions, genes)
