
import config

import sys

import numpy

from variation.variations import VariationsH5

import passport
import pop_building
import genome
import introgressions
import kegg


def _calc_gene_introgression_freqs(solgenomics_gene_id, genes, vars_index, introgression_freqs):
    gene = genes.get_gene(solgenomics_gene_id)
    idx0 = vars_index.index_pos(gene['Chromosome'].encode(), gene['Start'])
    idx1 = vars_index.index_pos(gene['Chromosome'].encode(), gene['End'])
    mask = numpy.logical_and(introgression_freqs.index >= idx0,
                             introgression_freqs.index <= idx1)
    gene_introgression_freqs = introgression_freqs.loc[mask].values
    return gene_introgression_freqs


def characterize_introgressions_per_pathways(variations, introgression_freqs, pathways, method, genes):

    vars_index = variations.pos_index

    mean_pathway_introgression_freqs = {}
    num_vars_with_most_introgressions = {}
    mean_gene_introgression_freqs = {}
    for pathway, info in pathways.items():

        gene_ids = info['solgenomic_gene_ids']
        pathway_introgression_freqs = []
        for gene_id in gene_ids:
            try:
                gene_introgression_freqs = _calc_gene_introgression_freqs(gene_id, genes, vars_index, introgression_freqs)
            except KeyError:
                continue

            pathway_introgression_freqs.extend(gene_introgression_freqs)
            mean_gene_introgression_freqs[gene_id] = numpy.mean(gene_introgression_freqs)

        if pathway_introgression_freqs:
            vars_with_most_introgression = numpy.array([])
            pathway_introgression_freqs = numpy.array(pathway_introgression_freqs)
            if method == 'quartile':
                mean_pathway_introgression_freqs[pathway] = numpy.nanpercentile(pathway_introgression_freqs, 90)
            elif method == 'mean':
                mean_pathway_introgression_freqs[pathway] = numpy.nanmean(pathway_introgression_freqs)
            elif method == 'mean_highest':
                q90 = numpy.nanpercentile(pathway_introgression_freqs, 90)
                vars_with_most_introgression = pathway_introgression_freqs[pathway_introgression_freqs > q90]
                num_vars_with_most_introgressions[pathway] = vars_with_most_introgression.size
                mean_pathway_introgression_freqs[pathway] = numpy.nanmean(vars_with_most_introgression)

    results = []
    pathways_sorted_by_introgression_freq = list(reversed(sorted(mean_pathway_introgression_freqs, key=lambda x: mean_pathway_introgression_freqs[x])))
    for pathway in pathways_sorted_by_introgression_freq:
        genes_in_pathway = list(pathways[pathway]['solgenomic_gene_ids'])
        num_genes_in_pathway = len(genes_in_pathway)
        result = {'kegg_path_id': pathway,
                  'mean_pathway_introgression_freqs': mean_pathway_introgression_freqs[pathway],
                  'num_vars_with_most_introgressions': num_vars_with_most_introgressions.get(pathway),
                  'num_genes_in_pathway': num_genes_in_pathway,
                  'description': pathways[pathway]['description']}

        genes_with_known_introgression_freq_in_pathway = [gene for gene in genes_in_pathway if mean_gene_introgression_freqs.get(gene) is not None and not numpy.isnan(mean_gene_introgression_freqs[gene])]
        genes_with_known_introgression_freq_in_pathway.sort(key=lambda x: mean_gene_introgression_freqs[x])

        genes = []
        for gene_id in reversed(genes_with_known_introgression_freq_in_pathway):
            gene = {'gene_id': gene_id,
                    'mean_gene_introgression_freq': mean_gene_introgression_freqs[gene_id]}
            genes.append(gene)
        result['genes'] = genes
        results.append(result)
    return results


def write_introgression_pathways_result(results, genes, fhand):
    for pathway_result in results:
        fhand.write(f'{pathway_result["kegg_path_id"]} {pathway_result["description"]}\n')
        fhand.write(f'Genes with introgressions:\n')
        for gene in pathway_result['genes']:
            introgression_freq = gene['mean_gene_introgression_freq']
            gene_id = gene["gene_id"]
            
            function = genes.get_annotated_function(gene_id)
            if function is None:
                function = ''

            gene_annotation = genes.gene_annotations[gene_id]

            if introgression_freq > 0:
                fhand.write(f'\t{gene["gene_id"]}\t {function}\tintrogression_freq: {introgression_freq}\n')

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

    res2 = kegg.get_pathways(add_solgenomic_ids=True, cache_dir=config.CACHE_DIR)
    pathways = res2['pathways']

    founder_pop = 'slc_ma'

    target_and_introgression_source = [('slc_pe', config.RANK1, 'sp_pe'),
                                       ('slc_ec', config.RANK1, 'sp_ec'),
                                       #('slc_pe_n_400m', config.RANK2, 'sp_pe'),
                                       #('slc_pe_n_1000m', config.RANK2, 'sp_pe'),
                                       #('slc_ec_n_600m', config.RANK2, 'sp_ec'),
                                       #('slc_ec_c_800m', config.RANK2, 'sp_ec')
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

        res = characterize_introgressions_per_pathways(variations, introgression_freqs, pathways, method, genes)

        path = out_dir / f'pathways_with_introgressions_for_{target_pop}.txt'
        write_introgression_pathways_result(res, genes, path.open('wt'))