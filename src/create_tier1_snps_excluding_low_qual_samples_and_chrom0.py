import config

from snp_filtering import get_low_quality_samples, filter_variations
from variation.variations import VariationsH5, VariationsArrays

from tomato_genome import get_pericentromic_regions


def take_some_snps(snps, n_snps, n_chunks=100, kept_fields=None):
    if n_snps > snps.num_variations:
        print('All SNPs taken: ', snps.num_variations)
        return snps

    out_snps = VariationsArrays()
    snps_per_chunk = int(snps.num_variations / n_chunks)
    snps_to_keep_per_chunk = round(n_snps / n_chunks)

    for chunk_idx in range(n_chunks - 1):
        start = chunk_idx * snps_per_chunk
        end = start + snps_to_keep_per_chunk + 1
        chunk = snps.get_chunk(slice(start, end), kept_fields=kept_fields)
        out_snps.put_chunks([chunk])
    print('Some SNPs taken: ', out_snps.num_variations)
    return out_snps


if __name__ == '__main__':
    
    debug = False
    use_low_qual_cache = False

    orig_vars = VariationsH5(str(config.ORIG_H5), 'r')
    out_variations = VariationsH5(str(config.TIER1_H5), 'w')

    filtering_stats_dir = config.FILTERING_STATS_DIR
    filtering_stats_dir.mkdir(exist_ok=True)
    tier1_stats_dir = filtering_stats_dir / 'tier1'
    tier1_stats_dir.mkdir(exist_ok=True)
    low_qual_samples_stats_dir = filtering_stats_dir / 'low_qual_samples'
    low_qual_samples_stats_dir.mkdir(exist_ok=True)
    missing_rate_hist_path = low_qual_samples_stats_dir / 'sample_called_rate.png'

    min_gt_dp = config.TIER1_MIN_GT_DP
    min_snp_qual = config.TIER1_MIN_SNP_QUAL
    max_mac = config.TIER1_MAX_MAC
    min_called = config.TIER1_MIN_CALLED
    max_het = config.TIER1_MAX_HET
    max_het_min_call_dp = config.TIER1_MAX_HET_MIN_CALL_DP
    
    low_quality_samples_txt = config.LOW_QUAL_SAMPLES_090
    kept_fields = config.RELEVANT_FIELDS
    min_called_rate_for_quality_samples = config.LOW_QUAL_MIN_SAMPLE_CALLED_RATE
    n_bins = config.LOW_QUAL_N_BINS

    print('input num_vars: ', orig_vars.num_variations)
    print('input num_samples: ', len(orig_vars.samples))
    if debug:
        orig_vars = take_some_snps(orig_vars, 20000)

    if use_low_qual_cache:
        low_quality_samples = [line.rstrip('\n') for line in low_quality_samples_txt.open()]
    else:
        low_quality_samples = get_low_quality_samples(take_some_snps(orig_vars, 1000000,
                                                      kept_fields=kept_fields),
                                                      min_called_rate_for_quality_samples,
                                                      regions_to_remove=None,
                                                      n_bins=n_bins,
                                                      missing_rate_hist_path=missing_rate_hist_path,
                                                      chunk_size=None,
                                                      max_chunks_to_process=None)

        low_quality_samples = list(low_quality_samples['low_quality_samples'])
    print("low quality samples: ", len(low_quality_samples), low_quality_samples)

    chrom0 = (b'SL2.50ch00', )
    regions_to_remove = [chrom0]

    filter_variations(orig_vars,
                      cache_dir=None,
                      chunk_size=600,
                      samples_to_remove=list(low_quality_samples),
                      min_gt_dp_setter=min_gt_dp,
                      filter_out_vars_with_non_major_allele_count_le=max_mac,
                      min_called=min_called,
                      max_het=max_het,
                      min_call_for_het=max_het_min_call_dp,
                      regions_to_remove=regions_to_remove,
                      sampling_rate=None,
                      out_variations=out_variations,
                      max_chunks_to_process=None,
                      fix_duplicated_alleles=True,
                      remove_non_variable_snvs=True,
                      kept_fields=kept_fields,
                      hist_path=tier1_stats_dir,
                      verbose=True)

    out_variations.samples = [sample.lower() for sample in out_variations.samples]
    print("output samples: ", len(out_variations.samples))
