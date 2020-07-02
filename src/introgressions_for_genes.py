
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

        gene_introgression_freqs = get_genes_with_most_introgressions(variations, introgression_freqs, genes, num_genes=10)
        print(gene_introgression_freqs)

'''

Genes with most introgresions in SLC EC
Solyc06g062400
    Chloroplast unusual positioning 1A
Solyc06g062900
    Transcription initiation factor TFIID subunit 12
Solyc06g062765
Solyc06g062760
    50S ribosomal protein L34
Solyc06g062590
    Transmembrane 9 superfamily protein member 4
Solyc06g062750
    Survival protein SurE-like phosphatase/nucleotidase 
Solyc06g062730
    Heterogeneous nuclear ribonucleoprotein A3 
Solyc09g074110
    AGAP009276
Solyc09g074360
    RNA binding protein
Solyc09g074330
    NDX1 homeobox protein
    
Genes with most introgresions in SLC PE
Solyc07g061780
    Ubiquitin carboxyl-terminal hydrolase family protein
Solyc03g005750
     Protein of unknown function
Solyc03g005680
    RNA-binding protein-like
Solyc11g072760
Solyc11g072750
    Chitinase a
    Like cellulose, chitin is an abundant biopolymer that is relatively resistant to degradation
Solyc11g072700
    Glycosyltransferase-like protein
Solyc11g072690
	kdsa
    3-desoxy-D-manno octulosonic acid-8-phosphate synthase
    The tomato kdsA gene expression and the relevant Kdo-8-P synthase activity were preferentially associated to dividing cells
Solyc07g062570
    Ubiquitin-conjugating enzyme E2 N
Solyc07g061920
    Glucan synthase like 3 
Solyc01g111560
    Ras-like GTP-binding protein RHO

'''
