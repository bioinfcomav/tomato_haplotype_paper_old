
import config

import pandas

from variation.variations import VariationsH5

import region
import passport
import pop_building
import genome


def analyze_regions(variations, regions_to_analyze, pops, pop_order, out_dir,
                    min_num_var_for_net, genes, remove_indels=False):

    for region_to_analyze in regions_to_analyze:
        reg_ = region_to_analyze
        this_region_out_dir = out_dir / f'{reg_["name"]}_{reg_["chrom"]}_{reg_["start"]}_{reg_["end"]}'
        this_region_out_dir.mkdir(exist_ok=True)
        print(this_region_out_dir)
        region.analyze_region(variations, region_to_analyze, pops, pop_order, this_region_out_dir,
                                min_num_var_for_net=min_num_var_for_net,
                                remove_indels=remove_indels, genes=genes, ignore_haplonet_failure=True)


def get_esther_regions():
    data = pandas.read_excel(OUT_DIR / 'genomic region for Jose.xlsx')

    regions = []
    for _, row in data.iterrows():
        chrom, start_end = row['region'].split(':')
        start_end = start_end.split('.')
        try:
            start, end = int(start_end[0]), int(start_end[-1])
        except ValueError:
            start_end = start_end[0].split('-')
            start, end = int(start_end[0]), int(start_end[-1])

        region = {'name': row['Name'],
                  'chrom': chrom,
                  'start': start, 'end': end}
        regions.append(region)
    return regions


if __name__ == '__main__':

    OUT_DIR = config.BASE_DIR /'esther'

    min_num_var_for_net = 14
    remove_indels = False

    variations = VariationsH5(str(config.TIER1_H5_LOWQ_085), 'r')
    passports = passport.get_sample_passports()
    pops = pop_building.get_pops({config.RANK1: config.ALL_POPS}, passports)

    regions_to_analyze = get_esther_regions()

    pop_order = ['sp_pe', 'sp_montane', 'sp_ec', 
                 'slc_ma',
                 'slc_pe', 'slc_ec',
                 'sll_mx',
                 'sll_vint', 'sll_modern']

    genes = genome.Genes()
    analyze_regions(variations, regions_to_analyze, pops, pop_order, OUT_DIR,
                    min_num_var_for_net=min_num_var_for_net, remove_indels=remove_indels,
                    genes=genes)
