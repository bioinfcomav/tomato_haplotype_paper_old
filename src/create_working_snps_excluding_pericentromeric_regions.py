import config

from snp_filtering import filter_variations
from variation.variations import VariationsH5
from solcap import determine_eucrohomatic_region


if __name__ == '__main__':
    
    input_vars = VariationsH5(str(config.TIER1_H5_LOWQ_085), 'r')
    output_vars = VariationsH5(str(config.WORKING_H5), 'w')

    min_called = config.TIER2_MIN_CALLED
    kept_fields = config.RELEVANT_FIELDS
    # pericentromeric_regions = [(chrom.encode(), limits[0], limits[1]) for chrom, limits in get_pericentromic_regions().items()]
    pericentromeric_regions =  [(chrom.encode(), start, end) for chrom, regions in euchromatic_regions.items() for start, end, is_euchromatic in regions if not is_euchromatic]

    filtering_stats_dir = config.FILTERING_STATS_DIR
    filtering_stats_dir.mkdir(exist_ok=True)
    working_set_stats_dir = filtering_stats_dir / 'working_set'
    working_set_stats_dir.mkdir(exist_ok=True)

    print('input num_vars: ', input_vars.num_variations)
    print('input num_samples: ', len(input_vars.samples))

    filter_variations(input_vars,
                      out_variations=output_vars,
                      cache_dir=None,
                      chunk_size=600,
                      min_called=min_called,
                      regions_to_remove=pericentromeric_regions,
                      remove_non_variable_snvs=True,
                      kept_fields=kept_fields,
                      hist_path=working_set_stats_dir,
                      verbose=True)

    print('output num_vars: ', output_vars.num_variations)
    print("output samples: ", len(output_vars.samples))
