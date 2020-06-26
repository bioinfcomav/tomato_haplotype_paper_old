
import config

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure


from variation.variations import VariationsH5

import genome
import passport
import pop_building
import allele_freq
import colors
import haplo_net


def get_variations_in_region(variations, chrom, start, end):

    chrom = chrom.encode()
    index = variations.pos_index
    idx0 = index.index_pos(chrom, start)
    idx1 = index.index_pos(chrom, end)
    return variations.get_chunk(slice(idx0, idx1))


if __name__ == '__main__':

    min_num_snp_for_network = 14
    remove_indels = False

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, passports)
    pop_for_samples = {sample: pop for pop, samples in pops.items() for sample in samples}

    out_dir = config.RELEVANT_GENES_DIR
    out_dir.mkdir(exist_ok=True)

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

    gene_infos = [edh2, rve8, prr95, prr73,
                  lhy, elf3, cab13, lux, elf4b, elf4a, gi,
                  prr37, toc1, cab2b, prr59, cab2, rve6, ccr2,
                  lnk1, phyb1, lnk2, fw32, fw113, fas, ovate, lc, fw22]

    genes = genome.Genes()

    for gene_info in gene_infos:
        print(gene_info['common_name'])

        gene = genes.get_gene(gene_info['solgenomics_id'])
        plus = gene_info.get('extra_nucleotides_around_for_freqs', 0)
        chrom = gene['Chromosome']
        start = gene['Start'] - plus
        end = gene['End'] + plus
        variations_for_freqs = get_variations_in_region(variations, chrom, start, end)
        print('n_vars_freqs', variations_for_freqs.num_variations)
        freqs = allele_freq.calc_allele_freq(variations_for_freqs, pops)

        fig = Figure()
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)

        allele_freq.plot_allele_freqs_bars(freqs, axes)

        plot_path = out_dir / f'{gene_info["common_name"]}.allele_freqs.svg'
        fig.tight_layout()
        fig.savefig(plot_path)


        plus = gene_info.get('extra_nucleotides_around_for_network', 0)
        chrom = gene['Chromosome']
        start = gene['Start'] - plus
        end = gene['End'] + plus
        variations_for_net = get_variations_in_region(variations, chrom, start, end)
        print('n_vars_net', variations_for_net.num_variations)

        res = haplo_net.calc_haplotype_network(variations_for_net,
                                               min_num_snp_for_network=min_num_snp_for_network,
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
