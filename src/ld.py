
import config

from pprint import pprint
import os
import hashlib, pickle, gzip
import shelve
import itertools

import numpy

from variation.variations import VariationsH5
from variation.variations.ld import calc_ld_along_genome

from snp_filtering import filter_variations
from passport import get_sample_passports
from pop_building import get_pops


def _yield_ld_cache(cache_path):
    with shelve.open(str(cache_path)) as ld_cache:
        for ld, dist, loci in ld_cache.values():
            yield ld, dist, loci


def calc_lds(variations, samples, max_dist, max_maf=0.95, chunk_size=100,
             min_num_samples=15, cache_dir=None):

    if len(samples) < min_num_samples:
        raise ValueError('Not enough samples')

    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += str(variations.num_variations)
        key += ','.join(sorted(samples))
        key += str(max_dist)
        key += str(max_maf)
        key += str(chunk_size)
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('lds_' + key + '.shelve')
        cache_db_fpath = str(cache_path) + '.db'
        if os.path.exists(cache_db_fpath):
            return _yield_ld_cache(cache_path)

    variations = filter_variations(variations,
                                   samples_to_keep=samples,
                                   cache_dir=cache_dir,
                                   max_maf=max_maf)
    print(f'calculating ld for {variations.num_variations} variations')
    ld_results = calc_ld_along_genome(variations,
                                      max_dist=max_dist,
                                      chunk_size=100,
                                      max_maf=max_maf)
    if cache_dir:
        with shelve.open(str(cache_path)) as ld_cache:
            for idx, (ld, dist, loci) in enumerate(ld_results):
                ld_cache[str(idx)] = (ld, dist, loci)
            ld_cache.close()
        return _yield_ld_cache(cache_path)
    else:
        return ld_results


def calc_lds_around_vars(variations, samples, win_size, max_maf=0.95,
                         min_num_lds_to_calc_ld=20,
                         max_dist_for_ld_calc=None,
                         min_num_samples=15,
                         cache_dir=None):
    if max_dist_for_ld_calc is None:
        max_dist_for_ld_calc = max_dist

    if win_size / 2 > max_dist_for_ld_calc:
        raise ValueError('win_size / 2 should be le than max_dist_for_ld_calc')

    if len(samples) < min_num_samples:
        raise ValueError('Not enough samples')

    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += str(variations.num_variations)
        key += ','.join(sorted(samples))
        key += str(win_size)
        key += str(max_maf)
        key += str(min_num_lds_to_calc_ld)
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('lds_for_vars' + key + '.pickle')
        if cache_path.exists():
            return pickle.load(gzip.open(cache_path, 'rb'))

    ld_results = calc_lds(variations, samples,
                          max_dist=max_dist_for_ld_calc,
                          max_maf=max_maf, min_num_samples=15,
                          cache_dir=cache_dir)

    half_win_size = win_size / 2

    vars_idx_by_pos = {}
    ld_sums_for_vars = numpy.zeros((variations.num_variations,))
    lds_seen_for_vars = numpy.zeros((variations.num_variations, ))
    for ld_result in ld_results:
        ld, dist, (chrom1, pos1, chrom2, pos2) = ld_result

        try:
            var1_idx = vars_idx_by_pos[(chrom1, pos1)]
        except KeyError:
            var1_idx = variations.pos_index.index_pos(chrom1, pos1)
            vars_idx_by_pos[(chrom1, pos1)] = var1_idx

        try:
            var2_idx = vars_idx_by_pos[(chrom2, pos2)]
        except KeyError:
            var2_idx = variations.pos_index.index_pos(chrom2, pos2)
            vars_idx_by_pos[(chrom2, pos2)] = var2_idx
        if chrom1 != chrom2:
            continue
        if abs(pos1 - pos2) > half_win_size:
            continue

        ld_sums_for_vars[var1_idx] += ld
        ld_sums_for_vars[var2_idx] += ld
        lds_seen_for_vars[var1_idx] += 1
        lds_seen_for_vars[var2_idx] += 1

    with numpy.errstate(invalid='ignore'):
        lds_for_vars = ld_sums_for_vars / lds_seen_for_vars

    lds_for_vars[lds_seen_for_vars < min_num_lds_to_calc_ld] = numpy.nan

    if cache_dir:
        pickle.dump(lds_for_vars, gzip.open(cache_path, 'wb'))

    return lds_for_vars


def _get_lds(variations, samples, max_maf,
             max_dist,
             min_num_samples,
             max_num_lds,
             cache_dir=None):

    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += str(variations.num_variations)
        key += ','.join(sorted(samples))
        key += str(max_dist)
        key += str(max_maf)
        key += str(max_num_lds)
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('ld_arrays_' + key + '.pickle.gz')
        if cache_path.exists():
            return pickle.load(gzip.open(cache_path, 'rb'))

    lds = calc_lds(variations, samples, max_maf=max_maf,
                   max_dist=max_dist,
                   min_num_samples=min_num_samples_in_pop,
                   cache_dir=cache_dir)
    lds = itertools.islice(lds, max_num_lds)
    lds, dists = zip(*((ld, dist) for ld, dist, _ in lds))
    lds = numpy.array(lds)
    dists = numpy.array(dists)

    mask = numpy.logical_not(numpy.isnan(dists))
    lds = lds[mask]
    dists = dists[mask]

    res = {'lds': lds, 'dists': dists}

    if cache_dir:
        pickle.dump(res, gzip.open(cache_path, 'wb'))

    return res


if __name__ == '__main__':

    max_maf = 0.95
    rank = config.RANK1
    max_dist = 1e5
    samples_for_loess = 10000
    num_lds_to_calculate = 50000
    min_num_samples_in_pop = 15
    cache_dir = None

    sample_passports = get_sample_passports()
    pops_descriptions = {rank: config.ALL_POPS}
    pops = get_pops(pops_descriptions, sample_passports)

    vars_path = config.WORKING_PHASED_H5
    variations = VariationsH5(str(vars_path), 'r')

    out_dir = config.LD_DIR / f'{rank[1]}_max_maf_{max_maf}_max_dist_{max_dist}'
    os.makedirs(out_dir, exist_ok=True)

    #variations = filter_variations(variations, max_maf=max_maf)

    pops = {pop:samples for pop, samples in pops.items() if len(samples) > min_num_samples_in_pop}

    for pop, samples in pops.items():
        lds = calc_lds(variations, samples, max_maf=max_maf,
                       max_dist=max_dist,
                       min_num_samples=min_num_samples_in_pop,
                       cache_dir=config.CACHE_DIR)
