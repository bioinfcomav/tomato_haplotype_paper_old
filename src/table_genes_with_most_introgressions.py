
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


def get_genes_with_most_introgressions_for_target_pops(variations, target_and_introgression_source, pops):

    genes_with_most_introgessions = {}
    for target_pop, target_pop_rank, introgression_source_pop in target_and_introgression_source:
        print(f'Introgressions in {target_pop}')
        samples_in_pop_with_possible_introgressions = pops[target_pop]
        samples_in_founder_pop = pops[founder_pop]
        samples_in_introgression_source_pop = pops[introgression_source_pop]

        introgession_config = {'samples_in_pop_with_introgressions': samples_in_pop_with_possible_introgressions,
                               'samples_in_founder_pop': samples_in_founder_pop,
                               'samples_in_introgression_source_pop': samples_in_introgression_source_pop,
                               'freq_threshold_to_consider_allele_present_in_founder_pop' : config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_PRESENT,
                               'freq_threshold_to_consider_allele_common_in_introgression_source_pop': config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_COMMON,
                               }
        introgression_freqs = introgressions.calc_introgession_freq_for_vars(variations, introgession_config)

        gene_introgression_freqs = get_genes_with_most_introgressions(variations,
                                                                      introgression_freqs,
                                                                      genes, num_genes=100)
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
    vars_path = config.WORKING_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, passports)

    win_size = 10000
    method = 'mean_highest'

    genes = genome.Genes()

    founder_pop = 'slc_ma'

    target_and_introgression_source = [('slc_pe_n', config.RANK1, 'sp_pe'),
                                       ('slc_pe_s', config.RANK1, 'sp_pe'),
                                       ('slc_ec', config.RANK1, 'sp_ec'),
                                       ]
    out_dir = config.FIGURES_DIR
    out_dir.mkdir(exist_ok=True)

    if True:
        genes_with_most_introgessions = get_genes_with_most_introgressions_for_target_pops(variations, target_and_introgression_source, pops)
    else:
        genes_with_most_introgessions = {'slc_pe_n': OrderedDict([('Solyc07g061780', 0.5609756097560976), ('Solyc03g005750', 0.55), ('Solyc03g005680', 0.55), ('Solyc11g072760', 0.5), ('Solyc11g072750', 0.5), ('Solyc11g072700', 0.4878048780487805), ('Solyc11g072690', 0.4878048780487805), ('Solyc07g062570', 0.4625), ('Solyc07g061920', 0.4605263157894737), ('Solyc01g111560', 0.45), ('Solyc09g065870', 0.43902439024390244), ('Solyc01g111050', 0.41025641025641024), ('Solyc07g062660', 0.40540540540540543), ('Solyc07g062040', 0.4), ('Solyc03g031900', 0.39473684210526316), ('Solyc03g031880', 0.39473684210526316), ('Solyc03g005060', 0.3902439024390244), ('Solyc03g005040', 0.3902439024390244), ('Solyc07g061980', 0.3875), ('Solyc07g062860', 0.38461538461538464), ('Solyc11g006470', 0.3780487804878049), ('Solyc11g006460', 0.3780487804878049), ('Solyc07g062160', 0.3780487804878049), ('Solyc07g062140', 0.3780487804878049), ('Solyc11g068760', 0.375), ('Solyc03g005670', 0.375), ('Solyc11g071650', 0.3717948717948718), ('Solyc04g079800', 0.3717948717948718), ('Solyc11g007710', 0.3684210526315789), ('Solyc11g007690', 0.3684210526315789), ('Solyc07g062940', 0.36585365853658536), ('Solyc07g062930', 0.36585365853658536), ('Solyc11g071730', 0.3625), ('Solyc11g070170', 0.3625), ('Solyc11g070160', 0.3625), ('Solyc10g007635', 0.358974358974359), ('Solyc02g071330', 0.35135135135135137), ('Solyc07g056660', 0.35), ('Solyc07g056650', 0.35), ('Solyc11g071700', 0.34615384615384615), ('Solyc10g007550', 0.34146341463414637), ('Solyc02g071250', 0.34146341463414637), ('Solyc02g071240', 0.34146341463414637), ('Solyc02g071140', 0.34146341463414637), ('Solyc11g008780', 0.3375), ('Solyc11g008180', 0.3333333333333333), ('Solyc10g006440', 0.3333333333333333), ('Solyc06g062400', 0.32926829268292684), ('Solyc08g066790', 0.32894736842105265), ('Solyc09g010950', 0.325), ('Solyc02g071340', 0.325), ('Solyc02g071320', 0.325), ('Solyc02g071260', 0.325), ('Solyc04g079840', 0.3170731707317073), ('Solyc08g066360', 0.3125), ('Solyc11g072430', 0.3108108108108108), ('Solyc11g072930', 0.3076923076923077), ('Solyc11g150145', 0.3076923076923077), ('Solyc11g068730', 0.3076923076923077), ('Solyc08g066030', 0.3055555555555556), ('Solyc08g065990', 0.3055555555555556), ('Solyc11g072820', 0.3048780487804878), ('Solyc08g066820', 0.3048780487804878), ('Solyc11g072730', 0.3), ('Solyc11g069160', 0.3), ('Solyc11g006620', 0.3), ('Solyc09g011030', 0.3), ('Solyc09g010990', 0.3), ('Solyc08g066110', 0.3), ('Solyc11g008720', 0.2948717948717949), ('Solyc11g008670', 0.2948717948717949), ('Solyc11g008610', 0.2948717948717949), ('Solyc10g007770', 0.2926829268292683), ('Solyc01g110470', 0.2894736842105263), ('Solyc11g007610', 0.2875), ('Solyc11g006490', 0.2875), ('Solyc10g006460', 0.2875), ('Solyc12g055920', 0.28205128205128205), ('Solyc12g055910', 0.28205128205128205), ('Solyc10g006030', 0.28205128205128205), ('Solyc02g093500', 0.28205128205128205), ('Solyc02g093460', 0.28205128205128205), ('Solyc02g070550', 0.28205128205128205), ('Solyc02g070520', 0.28205128205128205), ('Solyc08g066530', 0.27631578947368424), ('Solyc11g008905', 0.275), ('Solyc11g008900', 0.275), ('Solyc11g008520', 0.275), ('Solyc09g074270', 0.275), ('Solyc02g093510', 0.275), ('Solyc01g091460', 0.275), ('Solyc01g111620', 0.2702702702702703), ('Solyc01g110010', 0.2702702702702703), ('Solyc01g109990', 0.2702702702702703), ('Solyc02g092780', 0.2692307692307692), ('Solyc11g068990', 0.2682926829268293), ('Solyc11g068980', 0.2682926829268293), ('Solyc11g008440', 0.2682926829268293), ('Solyc01g110970', 0.2682926829268293), ('Solyc01g110460', 0.2682926829268293)]), 'slc_ec': OrderedDict([('Solyc06g062400', 0.9024390243902439), ('Solyc06g062900', 0.8780487804878049), ('Solyc06g062765', 0.8780487804878049), ('Solyc06g062760', 0.8780487804878049), ('Solyc06g062590', 0.8780487804878049), ('Solyc06g062750', 0.875), ('Solyc06g062730', 0.875), ('Solyc09g074110', 0.8536585365853658), ('Solyc09g074360', 0.8333333333333334), ('Solyc09g074330', 0.8333333333333334), ('Solyc06g063290', 0.8333333333333334), ('Solyc06g062610', 0.8333333333333334), ('Solyc06g062600', 0.8333333333333334), ('Solyc01g102370', 0.8333333333333334), ('Solyc06g063310', 0.8313008130081301), ('Solyc06g072760', 0.8292682926829268), ('Solyc06g063200', 0.8292682926829268), ('Solyc06g063190', 0.8292682926829268), ('Solyc06g082120', 0.825), ('Solyc06g082100', 0.825), ('Solyc01g104070', 0.825), ('Solyc01g104060', 0.825), ('Solyc06g082090', 0.8095238095238095), ('Solyc06g076965', 0.8095238095238095), ('Solyc01g104030', 0.8072289156626506), ('Solyc06g063050', 0.8055555555555556), ('Solyc01g090210', 0.782051282051282), ('Solyc12g014380', 0.7804878048780488), ('Solyc06g082200', 0.7804878048780488), ('Solyc06g065170', 0.7804878048780488), ('Solyc06g065060', 0.7804878048780488), ('Solyc01g090070', 0.7804878048780488), ('Solyc06g064830', 0.775), ('Solyc01g090200', 0.775), ('Solyc01g090170', 0.775), ('Solyc06g063250', 0.7692307692307693), ('Solyc06g063240', 0.7692307692307693), ('Solyc01g090230', 0.7692307692307693), ('Solyc06g065580', 0.7682926829268293), ('Solyc06g076670', 0.7625), ('Solyc06g076680', 0.7625), ('Solyc06g076650', 0.7625), ('Solyc01g102270', 0.7625), ('Solyc01g110460', 0.7619047619047619), ('Solyc01g111620', 0.7564102564102564), ('Solyc01g111590', 0.7564102564102564), ('Solyc12g015650', 0.7560975609756098), ('Solyc12g014615', 0.7560975609756098), ('Solyc08g067420', 0.7560975609756098), ('Solyc06g064660', 0.7560975609756098), ('Solyc06g063210', 0.7560975609756098), ('Solyc01g111190', 0.7560975609756098), ('Solyc01g110550', 0.7560975609756098), ('Solyc01g112160', 0.7532012195121951), ('Solyc06g084330', 0.75), ('Solyc06g084250', 0.75), ('Solyc02g069130', 0.7439024390243902), ('Solyc01g112350', 0.7439024390243902), ('Solyc06g084210', 0.7380952380952381), ('Solyc06g071900', 0.7380952380952381), ('Solyc06g071860', 0.7380952380952381), ('Solyc01g112240', 0.7375), ('Solyc01g112210', 0.7375), ('Solyc01g111200', 0.7317073170731707), ('Solyc01g110120', 0.7317073170731707), ('Solyc12g014540', 0.7307692307692307), ('Solyc12g014460', 0.7307692307692307), ('Solyc01g090100', 0.725), ('Solyc06g084460', 0.7195121951219512), ('Solyc06g084360', 0.7195121951219512), ('Solyc01g111150', 0.7195121951219512), ('Solyc01g111130', 0.7195121951219512), ('Solyc06g084080', 0.7142857142857143), ('Solyc06g083920', 0.7142857142857143), ('Solyc06g083870', 0.7142857142857143), ('Solyc06g071980', 0.7142857142857143), ('Solyc04g079800', 0.7125), ('Solyc01g111060', 0.7125), ('Solyc12g015660', 0.7073170731707317), ('Solyc06g071460', 0.7073170731707317), ('Solyc06g065240', 0.7073170731707317), ('Solyc06g065180', 0.7073170731707317), ('Solyc01g111760', 0.7051282051282052), ('Solyc01g111090', 0.6951219512195121), ('Solyc01g111050', 0.6951219512195121), ('Solyc09g090100', 0.6904761904761905), ('Solyc06g082980', 0.6904761904761905), ('Solyc02g062620', 0.6904761904761905), ('Solyc02g062570', 0.6904761904761905), ('Solyc01g110280', 0.6904761904761905), ('Solyc01g111560', 0.6875), ('Solyc01g111780', 0.6858974358974359), ('Solyc06g083550', 0.6829268292682927), ('Solyc06g083540', 0.6829268292682927), ('Solyc06g069510', 0.6829268292682927), ('Solyc04g077930', 0.6794871794871795), ('Solyc04g078040', 0.6785714285714286), ('Solyc09g072970', 0.675), ('Solyc08g068050', 0.675), ('Solyc02g068930', 0.6710526315789473)])}

    path = config.TABLE_GENES_WITH_MANY_INTROGRESSIONS
    write_genes_with_introgressions_table(path.open('wt'), genes_with_most_introgessions, genes)
