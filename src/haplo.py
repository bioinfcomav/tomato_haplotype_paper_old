
import numpy
import pandas

from variation import GT_FIELD, CHROM_FIELD, POS_FIELD, MISSING_INT
from variation.variations.filters import (SampleFilter, FLT_VARS,
                                          VariableAndNotAllMissing)


def get_haplos(variations, generate_uniq_haplos, generate_sample_haplos):
    samples = variations.samples

    gts = variations[GT_FIELD]
    ploidy = variations.ploidy

    haplotypes = None
    all_haplo_ids = []
    for ploidy_idx in range(ploidy):
        haplo_ids = [(sample, ploidy_idx) for sample in samples]
        all_haplo_ids.extend(haplo_ids)
        one_haploid_haplos = pandas.DataFrame(gts[:, :, ploidy_idx], columns=haplo_ids)
        if haplotypes is None:
            haplotypes = one_haploid_haplos
        else:
            haplotypes = pandas.concat([haplotypes, one_haploid_haplos], axis=1)

    res = {}

    if generate_sample_haplos:
        res['sample_haplos'] = haplotypes

    if generate_uniq_haplos:
        uniq_haplos, uniq_haplo_ids_that_correpond_to_all_haplos_idxs = numpy.unique(haplotypes.values, return_inverse=True, axis=1)
        uniq_haplos_ids = list(range(uniq_haplos.shape[1]))
        uniq_haplos = pandas.DataFrame(uniq_haplos, columns=uniq_haplos_ids)

        uniq_haplo_ids_that_correpond_to_all_haplos = pandas.Series(uniq_haplo_ids_that_correpond_to_all_haplos_idxs, index=all_haplo_ids)

        res['uniq_haplos'] = uniq_haplos
        res['uniq_haplo_ids_that_correpond_to_all_haplos'] = uniq_haplo_ids_that_correpond_to_all_haplos

    return res


def get_uniq_haplos(variations, generate_uniq_haplos, generate_sample_haplos):
    return get_haplos(variations, generate_uniq_haplos=True, generate_sample_haplos=False)


def _generate_haplos_along_genome(variations, win_params,
                                  generate_uniq_haplos,
                                  generate_sample_haplos,
                                  keep_only_variable_snps_in_samples,
                                  num_wins_to_process=None,
                                  samples=None):

    if samples is None:
        filter_samples = None
    else:
        sample_filter_object = SampleFilter(samples)
        filter_samples = lambda vars: sample_filter_object(vars)[FLT_VARS]

    if keep_only_variable_snps_in_samples:
        all_sample_filter_object = SampleFilter(keep_only_variable_snps_in_samples)
        variable_filter_object = VariableAndNotAllMissing()
        filter_non_variable = lambda vars: variable_filter_object(all_sample_filter_object(vars)[FLT_VARS])[FLT_VARS]
    else:
        filter_non_variable = None

    chunks = variations.iterate_wins(win_params['win_size'])

    n_processesed_chunks = 0
    for chunk_idx, chunk in enumerate(chunks):
        if chunk.num_variations < win_params['min_num_snp_for_window']:
                continue

        chrom = chunk[CHROM_FIELD][0].decode()
        chrom1 = chunk[CHROM_FIELD][-1].decode()
        if chrom != chrom1:
            continue
        start = chunk[POS_FIELD][0]
        #print(n_processesed_chunks, chunk_idx, chrom, start)

        if filter_non_variable:
            chunk = filter_non_variable(chunk)

            if chunk.num_variations < win_params['min_num_snp_for_window']:
                    continue

        if filter_samples:
            chunk = filter_samples(chunk)

        res = {'chrom': chrom, 'win_start': start}

        res.update(get_haplos(chunk, generate_uniq_haplos, generate_sample_haplos))

        yield res

        n_processesed_chunks += 1

        if num_wins_to_process and n_processesed_chunks >= num_wins_to_process:
            break


def generate_uniq_haplos_along_genome(variations, win_params, num_wins_to_process=None,
                                      keep_only_variable_snps_in_samples=None, samples=None):

    yield from _generate_haplos_along_genome(variations, win_params,
                                             generate_uniq_haplos=True,
                                             generate_sample_haplos=False,
                                             keep_only_variable_snps_in_samples=keep_only_variable_snps_in_samples,
                                             num_wins_to_process=num_wins_to_process,
                                             samples=samples)


def generate_sample_haplos_along_genome(variations, win_params, num_wins_to_process=None,
                                        keep_only_variable_snps_in_samples=None, samples=None):

    yield from _generate_haplos_along_genome(variations, win_params,
                                             generate_uniq_haplos=False,
                                             generate_sample_haplos=True,
                                             keep_only_variable_snps_in_samples=keep_only_variable_snps_in_samples,
                                             num_wins_to_process=num_wins_to_process,
                                             samples=samples)


def generate_haplos_along_genome(variations, win_params, num_wins_to_process=None,
                                 keep_only_variable_snps_in_samples=None, samples=None):

    yield from _generate_haplos_along_genome(variations, win_params,
                                             generate_uniq_haplos=True,
                                             generate_sample_haplos=True,
                                             keep_only_variable_snps_in_samples=keep_only_variable_snps_in_samples,
                                             num_wins_to_process=num_wins_to_process,
                                             samples=samples)


def generate_df_for_all_haplos_from_one_with_uniq_haplos_in_index(df_with_uniq_haplos_index, haplos_info):
    uniq_haplo_ids_that_correpond_to_all_haplos = haplos_info['uniq_haplo_ids_that_correpond_to_all_haplos']
    df_with_all_haplos_index = df_with_uniq_haplos_index.reindex(uniq_haplo_ids_that_correpond_to_all_haplos)
    df_with_all_haplos_index.index = uniq_haplo_ids_that_correpond_to_all_haplos.index

    df_with_all_haplos_index = df_with_all_haplos_index.dropna(how='all')

    return df_with_all_haplos_index


def filter_haplotypes_by_sample(haplos, haplos_ids, given_haplos_should_be_kept=False):
    if given_haplos_should_be_kept:
        return haplos.loc[:, haplos_ids]
    else:
        return haplos.drop(columns=haplos_ids)


def calc_edit_dist_between_indi_gts(gts_ind_pop1, gts_ind_pop2, win_size, min_num_snp_for_dist):
    debug = False

    if debug:
        print('gts_ind_pop1')
        print(gts_ind_pop1)
        print('gts_ind_pop2')
        print(gts_ind_pop2)
    missing_vars = numpy.logical_or(gts_ind_pop1 == MISSING_INT,
                                    gts_ind_pop2 == MISSING_INT)
    if gts_ind_pop1.shape[0] - numpy.sum(missing_vars) < min_num_snp_for_dist:
        return numpy.nan

    var_is_variable = numpy.logical_and(gts_ind_pop1 != gts_ind_pop2,
                                        numpy.logical_not(missing_vars))
    if debug:
        print(numpy.sum(var_is_variable), numpy.sum(var_is_variable) / win_size)

    dist = numpy.sum(var_is_variable) / win_size
    return dist


def calc_pairwise_dists_among_haplos(haplos, win_size, min_num_snp_for_dist):

    num_haplotypes = haplos.shape[1]

    dists = numpy.zeros((num_haplotypes, num_haplotypes))
    haplo_array = haplos.values
    for haplo_idx1 in range(haplo_array.shape[1]):
        haplo1 = haplo_array[:, haplo_idx1]
        for haplo_idx2 in range(haplo_idx1, haplo_array.shape[1]):
            haplo2 = haplo_array[:, haplo_idx2]

            dist = calc_edit_dist_between_indi_gts(haplo1,
                                                   haplo2,
                                                   win_size,
                                                   min_num_snp_for_dist)
            dists[haplo_idx1, haplo_idx2] = dist
            dists[haplo_idx2, haplo_idx1] = dist

    dists = pandas.DataFrame(dists, index=haplos.columns, columns=haplos.columns)
    return dists


def parse_haplo_id(haplo_id):
    chrom, win_start, sample, haploid_idx = haplo_id.split('%')
    win_start = int(win_start)
    haploid_idx = int(haploid_idx)
    return chrom, win_start, sample, haploid_idx


def create_haplo_id(chrom, win_start, sample, haploid_idx):
    return f'{chrom}%{win_start}%{sample}%{haploid_idx}'


def get_pop_classification_for_haplos(haplo_ids, pops):

    pops_for_samples = {sample: pop for pop, samples in pops.items() for sample in samples}
    pop_classification = {}
    for haplo_id in haplo_ids:
        sample = parse_haplo_id(haplo_id)[2]
        pop_classification[haplo_id] = pops_for_samples[sample]
    return pop_classification
