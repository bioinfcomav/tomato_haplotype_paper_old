
import config

from collections import defaultdict
from pprint import pprint
import time

import numpy
from pandas import DataFrame

from variation import CHROM_FIELD
from variation.variations import VariationsH5
from variation.variations.filters import VarsSamplingFilter2, FLT_VARS

from passport import get_sample_passports
from sample_selection import get_samples_for_criteria
from snp_filtering import (filter_variations,
                           keep_most_variable_snp_per_window,
                           keep_the_var_with_lowest_missing_gts_per_haplo_block)
import colors
from plot import plot_var_density_along_genome
from pop_distances import calc_kosman_dists
from pca import do_pcoa_from_dists, write_pca_curlywhirly_file


def get_sample_selection_criteria():

    rank1 = config.RANK1
    keep = config.KEEP
    remove = config.REMOVE

    criteria = []
    samples_to_remove = []
    samples_to_keep = ['hola']

    if True:
        criteria.append((rank1, ['sp_pe', 'sp_ec', 'sp_inter-andean'], KEEP))
        #criteria.append((RANK2, ['slc_ecu_big'], REMOVE))

    return {'criteria': criteria, 'samples_to_remove': samples_to_remove,
            'samples_to_keep': samples_to_keep}


def prepare_vars(variations,
                 plot_var_density=False,
                 win_size_to_keep_only_the_most_variable_var=None,
                 keep_at_most_n_vars=None,
                 difference_rate_allowed_for_haplo_block=0.2,
                 min_num_snps_to_use=100,
                 ignore_chromosome_representation_check=False,
                 cache_dir=None):

    if plot_var_density:
        plot_path = config.MULTIVAR_DIR / 'tier2_vars.svg'
        plot_var_density_along_genome(variations, plot_path=plot_path)

    if win_size_to_keep_only_the_most_variable_var:
        variations = keep_most_variable_snp_per_window(variations,
                                                       win_size=win_size_to_keep_only_the_most_variable_var,
                                                       cache_dir=cache_dir)
        print(f'one variable var per window: {variations.num_variations}')

        plot_path = config.MULTIVAR_DIR / 'tier2_vars_most_variable_per_region.svg'
        plot_var_density_along_genome(variations, plot_path=plot_path)

    if difference_rate_allowed_for_haplo_block:
        variations = keep_the_var_with_lowest_missing_gts_per_haplo_block(variations,
                                                                          difference_rate_allowed=difference_rate_allowed_for_haplo_block,
                                                                          cache_dir=cache_dir)

    if keep_at_most_n_vars:
        if variations.num_variations > keep_at_most_n_vars:
            variations = VarsSamplingFilter2(num_vars=keep_at_most_n_vars)(variations)[FLT_VARS]

        print(f'purged vars: {variations.num_variations}')

        plot_path = config.MULTIVAR_DIR / 'tier2_vars_purged.svg'
        plot_var_density_along_genome(variations, plot_path=plot_path)

    vars_for_dists = variations
    if vars_for_dists.num_variations > keep_at_most_n_vars:
        raise RuntimeError('Too many SNPs, some regions with high LD could be oversampled')
    if min_num_snps_to_use and vars_for_dists.num_variations  < min_num_snps_to_use:
        raise RuntimeError('Too few SNPs, maybe you should reduce the sampling rate')

    if numpy.unique(vars_for_dists[CHROM_FIELD]).size < 12:
        if ignore_chromosome_representation_check:
            warnings.warn(colors.TERMINAL_RED + 'Some chromosomes are not represented' + colors.TERMINAL_ENDC)
        else:
            raise RuntimeError('There are so few SNPs that some chromosomes are not represented')

    return vars_for_dists


def write_multivariant_result_for_curly(multivar_result, passports):
    field_paths_for_curly = [('classification', 'rank1'),
                        ('classification', 'rank2'),
                        ('country',),
                        #'morpho_type',
                        #'sw_group', 'tmp', 'region', 'het',
                        #'percent_haplos_close_to_ref',
                        #'sw_category'
                        ]

    passports_by_cat = defaultdict(dict)
    fields_for_curly = set()
    for sample_id, passport in passports.items():
        pprint(passport)
        for cat_path_in_passport in field_paths_for_curly:
            curly_cat = cat_path_in_passport[-1]

            passport_item = passport
            for key in cat_path_in_passport:
                print(passport_item)
                passport_item = passport_item.get(key, {})

            if passport_item:
                if isinstance(passport_item, dict):
                    print(passport_item)
                    raise ValueError('passport item should not be a dict, but a str')

                value = passport_item
                print('value', curly_cat, value)
                fields_for_curly.add(curly_cat)
                passports_by_cat[curly_cat][sample_id] = value

    passports_for_curly = {cat: samples for cat, samples in passports_by_cat.items() if cat in fields_for_curly}

    multivar_result['projections'] = DataFrame(multivar_result['projections'],
                                               index=multivar_result['samples'])

    multivar_dir = config.MULTIVAR_DIR
    multivar_dir.mkdir(exist_ok=True)
    curly_path =  multivar_dir / 'pcoa.curly'
    write_pca_curlywhirly_file(multivar_result,
                               curly_path,
                               categories=passports_for_curly)

    back_dir = multivar_dir / 'back'
    back_dir.mkdir(exist_ok=True)
    datestamp = time.strftime("%Y-%m-%d_%H:%M:%S", time.gmtime())
    curly_path_back = back_dir / f'pcoa.{datestamp}.pcoa'
    write_pca_curlywhirly_file(multivar_result,
                               curly_path_back,
                               categories=passports_for_curly)


if __name__ == '__main__':

    config.MULTIVAR_DIR.mkdir(exist_ok=True)

    debug = False

    filter_by_maf = True
    chunk_size = 1000

    plot_var_density = False
    win_size_to_keep_only_the_most_variable_var = 1e5
    keep_at_most_n_vars = 2000
    difference_rate_allowed_for_haplo_block = 0.2
    max_num_vars_for_pca = 2000
    min_num_snps_to_use = 100
    ignore_chromosome_representation_check = False

    if debug:
        max_chunks_to_process = 1000
        cache_dir = None
    else:
        max_chunks_to_process = None
        cache_dir=config.CACHE_DIR

    passports = get_sample_passports()

    vars_path = config.TIER1_PHASED_LOW_QUAL_09_MISSING_085
    variations = VariationsH5(str(vars_path), 'r')

    all_samples = variations.samples
    print(variations.num_variations)

    criteria = get_sample_selection_criteria()
    samples_to_use = get_samples_for_criteria(all_samples,
                                              passports,
                                              criteria,
                                              skip_samples_with_no_passport=config.SKIP_SAMPLES_WITH_NO_PASSPORT)
    
    if filter_by_maf:
        print(colors.TERMINAL_BLUE + 'Doing a MAF filtering you are reducing the distance of the samples that do not belong to a well represented population' + colors.TERMINAL_ENDC)
        max_maf = config.TIER2_PCA_MAX_MAF
    else:
        max_maf = None

    print(variations.num_variations)
    tier2_vars = filter_variations(variations,
                                   chunk_size=chunk_size,
                                   samples_to_keep=samples_to_use,
                                   cache_dir=cache_dir,
                                   max_mac=config.TIER2_MAX_MAC,
                                   max_maf=max_maf,
                                   min_called=config.TIER2_MIN_CALLED,
                                   max_het=config.TIER2_MAX_HET,
                                   min_call_for_het=config.TIER2_MAX_HET_MIN_CALL_DP,
                                   kept_fields=config.RELEVANT_FIELDS,
                                   max_chunks_to_process=max_chunks_to_process,
                                   remove_non_variable_snvs=True,
                                   verbose=True
                                   )

    if not tier2_vars.num_variations:
        raise ValueError('No SNPs left after tier2')

    print(f'tier2 vars: {tier2_vars.num_variations}')

    var_counts_per_chrom = dict(zip(*numpy.unique(tier2_vars[CHROM_FIELD], return_counts=True)))

    if not debug and len(var_counts_per_chrom) < 12:
        print('Only some chromosomes have SNPs:')
        for chrom, num_vars in var_counts_per_chrom.items():
            chrom = chrom.decode()
            print(f'\t{chrom}:\t{num_vars}')
        raise RuntimeError('Some chromosomes have no SNPs')

    vars_for_dists = prepare_vars(tier2_vars,
                                  plot_var_density=plot_var_density,
                                  win_size_to_keep_only_the_most_variable_var=win_size_to_keep_only_the_most_variable_var,
                                  keep_at_most_n_vars=keep_at_most_n_vars,
                                  difference_rate_allowed_for_haplo_block=difference_rate_allowed_for_haplo_block,
                                  min_num_snps_to_use=min_num_snps_to_use,
                                  ignore_chromosome_representation_check=ignore_chromosome_representation_check,
                                  cache_dir=cache_dir)
    print(f'Num vars for dists: {vars_for_dists.num_variations}')

    dist_result = calc_kosman_dists(vars_for_dists, cache_dir=config.CACHE_DIR)
    multivar_result = do_pcoa_from_dists(dist_result)

    if 'var_percentages' in multivar_result:
        tot_variance = sum(multivar_result['var_percentages'])
        variance = sum(multivar_result['var_percentages'][:3]) / tot_variance * 100
        print('Percentage of variance first 3 components: ',
            str(multivar_result['var_percentages'][:3]), ' (%.1f%%)' % variance)
    print('Samples in PCA: ', len(multivar_result['samples']))

    write_multivariant_result_for_curly(multivar_result, passports)
