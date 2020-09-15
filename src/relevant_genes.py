
import config

import math

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5
from variation.variations.filters import SNPPositionFilter, FLT_VARS

import genome
import passport
import pop_building
import imputation
import allele_freq
import haplo_net
import colors




if __name__ == '__main__':

    remove_indels = False

    min_num_vars_for_imputation = 100
    min_num_var_for_haplo_freqs = 5
    min_num_var_for_net = 14

    num_bp_to_increase_in_region = 0

    vars_path = config.TIER1_H5_LOWQ_085
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, passports)
    pop_for_samples = {sample: pop for pop, samples in pops.items() for sample in samples}

    pop_order = ['sp_pe', 'sp_montane', 'sp_ec', 
                 'slc_ma',
                 'slc_pe', 'slc_ec',
                 'sll_mx',
                 'sll_vint']

    out_dir = config.RELEVANT_GENES_DIR
    out_dir.mkdir(exist_ok=True)

    # rice genes related to latitude
    # https://www.frontiersin.org/articles/10.3389/fpls.2020.00864/full
    # Hd1, DTH8, Ghd7, OsCOL4, DTH7, Hd6, Se5, and PhyB
    # rice Hd1: OsBBX26
    # rice DTH8: Os08g0174500

    se5 = {'common_name': 'se5',
            'solgenomics_id': 'Solyc12g009470'} 

    phyba = {'common_name': 'phyba',
            'solgenomics_id': 'Solyc01g059870'} 
    phybb = {'common_name': 'phybb',
            'solgenomics_id': 'Solyc05g053410'} 
    phybc = {'common_name': 'phybc',
            'solgenomics_id': 'Solyc02g071260'} 
    phybd = {'common_name': 'phybd',
            'solgenomics_id': 'Solyc07g045480'} 
    phybe = {'common_name': 'phybe',
            'solgenomics_id': 'Solyc10g044670'} 

    hd6a = {'common_name': 'hd6a',
            'solgenomics_id': 'Solyc03g044140'} 
    hd6b = {'common_name': 'hd6b',
            'solgenomics_id': 'Solyc02g064700'} 
    hd6c = {'common_name': 'hd6c',
            'solgenomics_id': 'Solyc02g092640'} 

    hd1 = {'common_name': 'hd1',
            'solgenomics_id': 'Solyc05g024010'}

    dth8a = {'common_name': 'dth8a',
            'solgenomics_id': 'Solyc04g054150'}

    dth8b = {'common_name': 'dth8b',
            'solgenomics_id': 'Solyc12g006120'}

    dth8c = {'common_name': 'dth8c',
            'solgenomics_id': 'Solyc07g065500'}

    dth8d = {'common_name': 'dth8d',
            'solgenomics_id': 'Solyc12g027650'}

    dth8e = {'common_name': 'dth8e',
            'solgenomics_id': 'Solyc09g007290'}

    dth8f = {'common_name': 'dth8f',
            'solgenomics_id': 'Solyc01g067130'}

    dth8g = {'common_name': 'dth8g',
            'solgenomics_id': 'Solyc04g009520'}

    dth8h = {'common_name': 'dth8h',
            'solgenomics_id': 'Solyc09g074760'}

    dth8i = {'common_name': 'dth8i',
            'solgenomics_id': 'Solyc06g009010'}

    dth8j = {'common_name': 'dth8j',
            'solgenomics_id': 'Solyc04g049910'}

    dth8k = {'common_name': 'dth8k',
            'solgenomics_id': 'Solyc04g015060'}

    dth8l = {'common_name': 'dth8l',
            'solgenomics_id': 'Solyc05g005390'}

    dth8m = {'common_name': 'dth8m',
            'solgenomics_id': 'Solyc05g005380'}
    dth8n = {'common_name': 'dth8n',
            'solgenomics_id': 'Solyc10g009440'}

    dth8o = {'common_name': 'dth8o',
            'solgenomics_id': 'Solyc02g032185'}

    dth8p = {'common_name': 'dth8p',
            'solgenomics_id': 'Solyc01g099320'}
    dth8q = {'common_name': 'dth8q',
            'solgenomics_id': 'Solyc05g015550'}
    dth8q = {'common_name': 'dth8q',
            'solgenomics_id': 'Solyc05g005440'}
    dth8r = {'common_name': 'dth8r',
            'solgenomics_id': 'Solyc05g005360'}
    dth8s = {'common_name': 'dth8s',
            'solgenomics_id': 'Solyc05g005350'}
    dth8t = {'common_name': 'dth8t',
            'solgenomics_id': 'Solyc07g065580'}

    dth8u = {'common_name': 'dth8u',
            'solgenomics_id': 'Solyc07g065570'}
    dth8v = {'common_name': 'dth8v',
            'solgenomics_id': 'Solyc07g065567'}
    dth8w = {'common_name': 'dth8w',
            'solgenomics_id': 'Solyc07g065563'}
    dth8x = {'common_name': 'dth8x',
            'solgenomics_id': 'Solyc05g005370'}
    dth8y = {'common_name': 'dth8y',
            'solgenomics_id': 'Solyc02g032190'}
    dth8z = {'common_name': 'dth8z',
            'solgenomics_id': 'Solyc06g069310'}

    oscol4a = {'common_name': 'oscol4a',
            'solgenomics_id': 'Solyc12g096500'}
    oscol4b = {'common_name': 'oscol4b',
            'solgenomics_id': 'Solyc07g006630'}
    oscol4c = {'common_name': 'oscol4c',
            'solgenomics_id': 'Solyc02g089540'}
    oscol4d = {'common_name': 'oscol4d',
            'solgenomics_id': 'Solyc08g006530'}
    oscol4e = {'common_name': 'oscol4e',
            'solgenomics_id': 'Solyc02g089520'}


    fw22 = {'common_name': 'fw22',
            'solgenomics_id': 'Solyc02g090730',
            'extra_nucleotides_around_for_freqs': 50000,
            'extra_nucleotides_around_for_network': 150000}

    lc = {'common_name': 'lc', # wuschel
            'solgenomics_id': 'Solyc02g083950',
            'extra_nucleotides_around_for_freqs': 40000,
            'extra_nucleotides_around_for_network': 200000}

    ovate = {'common_name': 'ovate',
            'solgenomics_id': 'Solyc02g085500',
            'extra_nucleotides_around_for_freqs': 30000,
            'extra_nucleotides_around_for_network': 100000}

    fas = {'common_name': 'fas', # YABBY
          'solgenomics_id': 'Solyc11g071810',
          'extra_nucleotides_around_for_freqs': 15000,
          'extra_nucleotides_around_for_network': 60000}

    fw113 = {'common_name': 'fw113',
            'solgenomics_id': 'Solyc11g071940',
            'extra_nucleotides_around_for_freqs': 30000,
            'extra_nucleotides_around_for_network': 50000}

    fw32 = {'common_name': 'fw32',    #KLUH
            'solgenomics_id': 'Solyc03g114940',
            'extra_nucleotides_around_for_freqs': 0,
            'extra_nucleotides_around_for_network': 0}

    # Mutations in EID1 and LNK2 caused light-conditional clock deceleration during tomato domestication
    lnk2 = {'common_name': 'lnk2',
            'solgenomics_id': 'Solyc01g068560',
            'extra_nucleotides_around_for_freqs': 30000,
            'extra_nucleotides_around_for_network': 80000}

    # Mutations in EID1 and LNK2 caused light-conditional clock deceleration during tomato domestication
    phyb1 = {'common_name': 'phyb1',    # EID1
            'solgenomics_id': 'Solyc01g059870',
            'extra_nucleotides_around_for_freqs': 5900000,
            'extra_nucleotides_around_for_network': 5950000}

    lnk1 = {'common_name': 'lnk1',
            'solgenomics_id': 'Solyc01g105120',
            'extra_nucleotides_around_for_freqs': 7500,
            'extra_nucleotides_around_for_network': 80000}

    ccr2 = {'common_name': 'ccr2',
            'solgenomics_id': 'Solyc01g109660',
            'extra_nucleotides_around_for_freqs': 40000,
            'extra_nucleotides_around_for_network': 80000}

    rve6 = {'common_name': 'rve6',
            'solgenomics_id': 'Solyc01g095030',
            'extra_nucleotides_around_for_freqs': 10000,
            'extra_nucleotides_around_for_network': 80000}

    cab2 = {'common_name': 'cab2',
            'solgenomics_id': 'Solyc02g071010',
            'extra_nucleotides_around_for_freqs': 5000,
            'extra_nucleotides_around_for_network': 80000}

    prr59 = {'common_name': 'prr59',
            'solgenomics_id': 'Solyc03g081240',
            'extra_nucleotides_around_for_freqs': 25000,
            'extra_nucleotides_around_for_network': 80000}

    cab2b = {'common_name': 'cab2b',
            'solgenomics_id': 'Solyc03g005770',
            'extra_nucleotides_around_for_freqs': 5000,
            'extra_nucleotides_around_for_network': 80000}

    toc1 = {'common_name': 'toc1',
            'solgenomics_id': 'Solyc03g115770',
            'extra_nucleotides_around_for_freqs': 25000,
            'extra_nucleotides_around_for_network': 90000}

    prr37 = {'common_name': 'prr37',
            'solgenomics_id': 'Solyc04g049670',
            'extra_nucleotides_around_for_freqs': 9029500,
            'extra_nucleotides_around_for_network': 9200000}

    gi = {'common_name': 'gi',
            'solgenomics_id': 'Solyc04g071990',
            'extra_nucleotides_around_for_freqs': 0,
            'extra_nucleotides_around_for_network': 200000}

    elf4a = {'common_name': 'elf4a',
            'solgenomics_id': 'Solyc06g051660',
            'extra_nucleotides_around_for_freqs': 10000,
            'extra_nucleotides_around_for_network': 30000}

    elf4b = {'common_name': 'elf4b',
            'solgenomics_id': 'Solyc06g051680',
            'extra_nucleotides_around_for_freqs': 15000,
            'extra_nucleotides_around_for_network': 40000}

    lux = {'common_name': 'luxa',
            'solgenomics_id': 'Solyc06g005680',
            'extra_nucleotides_around_for_freqs': 3000,
            'extra_nucleotides_around_for_network': 30000}

    cab13 = {'common_name': 'cab13',
            'solgenomics_id': 'Solyc07g063600',
            'extra_nucleotides_around_for_freqs': 30000,
            'extra_nucleotides_around_for_network': 50000}

    elf3 = {'common_name': 'elf3',
            'solgenomics_id': 'Solyc08g065870',
            'extra_nucleotides_around_for_freqs': 30000,
            'extra_nucleotides_around_for_network': 80000}

    lhy = {'common_name': 'lhy',
            'solgenomics_id': 'Solyc10g005080',
            'extra_nucleotides_around_for_freqs': 1000,
            'extra_nucleotides_around_for_network': 200000}

    prr73 = {'common_name': 'prr73',
            'solgenomics_id': 'Solyc10g086000',
            'extra_nucleotides_around_for_freqs': 1000,
            'extra_nucleotides_around_for_network': 200000}

    prr95 = {'common_name': 'prr95',
            'solgenomics_id': 'Solyc10g005030',
            'extra_nucleotides_around_for_freqs': 10000,
            'extra_nucleotides_around_for_network': 200000}

    rve8 = {'common_name': 'rve8',
            'solgenomics_id': 'Solyc10g084370',
            'extra_nucleotides_around_for_freqs': 20000,
            'extra_nucleotides_around_for_network': 150000}

    edh2 = {'common_name': 'edh2',
            'solgenomics_id': 'Solyc09g075090',
            'extra_nucleotides_around_for_freqs': 3000,
            'extra_nucleotides_around_for_network': 100000}

    # phase QTL in Müller 2016
    # This polymorphism (a 3-bp deletion) is located in a region highly conserved among
    # several Solanaceae species and probably represents the causal mutation underlying the phase QTL
    # NDGNRKkIDLNCAFCSSK
    #  NDGNRKIDLNCAFCSSK
    edh1 = {'common_name': 'edh1',
            'solgenomics_id': 'Solyc09g075080'}
    edh1_mut = {'common_name': 'edh1_mut',
                'chrom': 'SL4.0ch09',
                'start': 63111577,
                'end': 63111579,
                'min_num_var_for_haplo_freqs': 2}

    gene_infos = [edh1, edh1_mut, edh2, rve8, prr95, prr73,
                  lhy, elf3, cab13, lux, elf4b, elf4a, gi,
                  prr37, toc1, cab2b, prr59, cab2, rve6, ccr2,
                  lnk1, phyb1, lnk2, fw32, fw113, fas, ovate, lc, fw22]

    genes = genome.Genes()

    #min_num_vars_for_imputation = 100
    #min_num_var_for_haplo_freqs = 5

    for gene_info in gene_infos:

        solgenomics_id = gene_info.get('solgenomics_id')

        common_gene_name = gene_info.get('common_name', solgenomics_id)
        print(common_gene_name)

        if solgenomics_id:
                gene = genes.get_gene(gene_info['solgenomics_id'])
                plus = gene_info.get('extra_nucleotides_around_for_freqs', 0)
                chrom = gene['Chromosome']
                gene_start = gene['Start'] - num_bp_to_increase_in_region
                if gene_start < 0:
                    gene_start = 0
                gene_end = gene['End'] + num_bp_to_increase_in_region
        else:
                chrom = gene_info['chrom']
                gene_start = gene_info['start']
                gene_end = gene_info['end']

        variations_in_region = get_variations_in_region(variations,
                                                        chrom, gene_start, gene_end,
                                                        min_num_vars_for_imputation)
        imputed_variations = imputation.impute_variations(variations_in_region,
                                                          genetic_map_path=config.BEAGLE_MAP)
        this_min_num_var_for_haplo_freqs = gene_info.get('min_num_var_for_haplo_freqs', min_num_var_for_haplo_freqs)
        variations_for_freqs = get_variations_in_region(imputed_variations, chrom, gene_start, gene_end,
                                                        this_min_num_var_for_haplo_freqs)
        print('n_vars_freqs', variations_for_freqs.num_variations)
        freqs = allele_freq.calc_allele_freq(variations_for_freqs, pops)

        fig = Figure()
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)

        allele_freq.plot_allele_freqs_bars(freqs, axes, pop_order=pop_order)

        plot_path = out_dir / f'{common_gene_name}.allele_freqs.svg'
        fig.tight_layout()
        fig.savefig(plot_path)

        variations_for_net = get_variations_in_region(imputed_variations, chrom, gene_start, gene_end,
                                                      min_num_var_for_net)
        
        print('n_vars_net', variations_for_net.num_variations)

        res = haplo_net.calc_haplotype_network(variations_for_net,
                                               min_num_snp_for_network=min_num_var_for_net,
                                               remove_indels=remove_indels)

        fig = Figure()
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)
        haplo_net.plot_network(res['network_nx'], res['sample_sets'],
                               axes, pop_for_samples,
                               pop_color_schema=colors.ColorSchema(colors.CLASSIFICATION_RANK1_COLORS))

        plot_path = out_dir / f'{gene_info["common_name"]}.network.svg'
        fig.tight_layout()
        fig.savefig(plot_path)
