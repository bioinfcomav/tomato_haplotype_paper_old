
import pandas

from variation.variations.multivariate import do_pcoa

from haplo import (generate_uniq_haplos_along_genome, filter_haplotypes_by_sample,
                   calc_pairwise_dists_among_haplos,
                   generate_df_for_all_haplos_from_one_with_uniq_haplos_in_index)


def do_pcoa_from_dists(dists):

    multivar_result = do_pcoa(dists.values)
    multivar_result['projections'] = pandas.DataFrame(multivar_result['projections'],
                                                      index=dists.index)

    return multivar_result


def _do_pcoa_for_haplos(haplos, win_size, min_num_snp_for_dist,
                        n_dims_to_keep, chrom, win_start, cache_dir=None):

    if cache_dir is not None:
        key = ''
        key += f'haplos_index_' + '-'.join(map(str,haplos.index))
        key += f'haplos_columns_' + '-'.join(map(str, haplos.columns))
        with numpy.printoptions(threshold=numpy.inf):
            key += str(haplos.values)
        key += f'win_size_{win_size}'
        key += f'min_num_snp_for_window_{min_num_snp_for_dist}'
        key += f'n_dims_to_keep_{n_dims_to_keep}'

        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('haplo_pcoa_' + key + '.pickle')
        if cache_path.exists():
            return pickle.load(cache_path.open('rb'))

    dists = calc_pairwise_dists_among_haplos(haplos, win_size=win_size, min_num_snp_for_dist=min_num_snp_for_dist)

    pcoa = do_pcoa_from_dists(dists)

    if n_dims_to_keep:
        pcoa['projections'] = pcoa['projections'].iloc[:, :n_dims_to_keep]

    if cache_dir is not None:
        pickle.dump(pcoa, cache_path.open('wb'))

    return pcoa


def do_pcoas_along_the_genome(variations, win_params, num_wins_to_process=None,
                              haplotypes_to_exclude=None,
                              haplotypes_to_include=None,
                              samples=None, n_dims_to_keep=None, cache_dir=None):

    if haplotypes_to_exclude is None:
        haplotypes_to_exclude = {}

    if haplotypes_to_include is None:
        haplotypes_to_include = {}

    pcoas = []
    for haplos_info in generate_uniq_haplos_along_genome(variations, win_params,
                                                         num_wins_to_process=num_wins_to_process,
                                                         samples=samples):
        haplos = haplos_info['uniq_haplos']
        win_size = win_params['win_size']
        min_num_snp_for_dist = win_params['min_num_snp_for_window']
        chrom = haplos_info['chrom']
        win_start = haplos_info['win_start']

        uniq_haplo_ids_that_correpond_to_all_haplos = haplos_info['uniq_haplo_ids_that_correpond_to_all_haplos']

        haplos_to_include_in_this_win = haplotypes_to_include.get((chrom, win_start))
        if not haplos_to_include_in_this_win and haplotypes_to_include:
            continue
        if haplos_to_include_in_this_win:
            uniq_haplos_to_include = {uniq_haplo_ids_that_correpond_to_all_haplos[(sample, haploid_idx)] for sample, haploid_idx in haplos_to_include_in_this_win}
            haplos = filter_haplotypes_by_sample(haplos, uniq_haplos_to_include, given_haplos_should_be_included=True)

        haplos_to_exclude_in_this_win = haplotypes_to_exclude.get((chrom, win_start))
        if haplos_to_exclude_in_this_win:
            uniq_haplos_to_exclude = {uniq_haplo_ids_that_correpond_to_all_haplos[(sample, haploid_idx)] for sample, haploid_idx in haplos_to_exclude_in_this_win}
            uniq_haplos_to_exclude = uniq_haplos_to_exclude.intersection(haplos.columns)
            haplos = filter_haplotypes_by_sample(haplos, uniq_haplos_to_exclude)

        if len(haplos.columns) < n_dims_to_keep:
            continue

        pcoa = _do_pcoa_for_haplos(haplos, win_size=win_size,
                                   min_num_snp_for_dist=min_num_snp_for_dist,
                                   n_dims_to_keep=n_dims_to_keep,
                                   chrom=chrom, win_start=win_start, cache_dir=cache_dir)

        projections_for_uniq_haplos = pcoa['projections']

        projections_for_all_haplos = generate_df_for_all_haplos_from_one_with_uniq_haplos_in_index(projections_for_uniq_haplos, haplos_info)

        if haplos_to_exclude_in_this_win:
            projections_for_all_haplos = projections_for_all_haplos.dropna(axis=0, how='any')

        pcoas.append({'projections': projections_for_all_haplos,
                      'chrom': haplos_info['chrom'],
                      'win_start': haplos_info['win_start']})

    return pcoas
