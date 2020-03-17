

import hashlib
import pickle
import random
import os

import numpy

from variation import GT_FIELD, CHROM_FIELD, POS_FIELD
from variation.variations.pipeline import Pipeline
from variation.variations.filters import (SampleFilter, FLT_VARS, COUNTS,
                                          EDGES, SNPPositionFilter,
                                          MacFilter, MinCalledGTsFilter,
                                          ObsHetFilter, VarsSamplingFilter,
                                          filter_samples_by_missing_rate,
                                          VariableAndNotAllMissing,
                                          MafFilter,
                                          LowDPGTsToMissingSetter)
from variation.plot import plot_histogram
from variation.variations import VariationsH5, VariationsArrays
from variation.variations.blocks import generate_blocks as variation_generate_blocks
from variation.variations.stats import calc_called_gt, calc_maf, calc_missing_gt


def render_histogram(res, fhand, vlines, hrange=None):
    xlower, xupper = min(res[EDGES]), max(res[EDGES]) if hrange == None else hrange
    plot_histogram(res[COUNTS], res[EDGES],
                    fhand=fhand,
                    vlines=[vlines],
                    mpl_params={'set_yscale': {'args': ['log']},
                                'set_xbound': {'kwargs': {'lower': xlower,
                                                            'upper': xupper}}})


def filter_variations(variations, samples_to_keep=None,
                      samples_to_remove=None, cache_dir=None,
                      min_gt_dp_setter=None,
                      chunk_size=600,
                      filter_out_vars_with_non_major_allele_count_le=None, 
                      max_maf=None, min_called=None,
                      max_het=None, regions_to_keep=None, min_call_for_het=0,
                      regions_to_remove=None, sampling_rate=None,
                      out_variations=None, max_chunks_to_process=None,
                      kept_fields=None, remove_non_variable_snvs=False,
                      hist_path=None, verbose=False,
                      ignore_sampling_rate_for_cache=False):

    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += str(variations.num_variations)
        key += str(samples_to_keep)
        key += str(samples_to_remove)
        key += str(max_mac)
        key += str(max_maf)
        key += str(min_called)
        key += str(max_het)
        key += str(regions_to_keep)
        key += str(regions_to_remove)
        key += str(min_call_for_het)
        if not ignore_sampling_rate_for_cache:
            key += str(sampling_rate)
        key = hashlib.md5(key.encode()).hexdigest()
        h5_cache_path = cache_dir / ('tier2_vars' + key + '.h5')
        #print(h5_cache_path)
        if h5_cache_path.exists():
            return VariationsH5(str(h5_cache_path), 'r')

    pipeline = Pipeline()

    if sampling_rate is not None:
        flt = VarsSamplingFilter(sampling_rate)
        pipeline.append(flt, id_='sampling')

    if samples_to_keep is not None:
        flt = SampleFilter(samples_to_keep)
        pipeline.append(flt, id_='samples_to_keep')

    if samples_to_remove is not None:
        flt = SampleFilter(samples_to_remove, reverse=True)
        pipeline.append(flt, id_='samples_to_remove')

    if regions_to_keep is not None:
        flt = SNPPositionFilter(regions_to_keep, reverse=False)
        pipeline.append(flt, id_='regions_to_keep')

    if regions_to_remove is not None:
        flt = SNPPositionFilter(regions_to_remove, reverse=True)
        pipeline.append(flt, id_='regions_to_remove')

    if min_gt_dp_setter is not None:
        flt = LowDPGTsToMissingSetter(min_dp=min_gt_dp_setter)
        pipeline.append(flt, id_='min_gt_dp_setter')

    if remove_non_variable_snvs:
        flt = VariableAndNotAllMissing()
        pipeline.append(flt, id_='variable_and_not_all_missing')

    do_histogram = bool(hist_path)

    if max_mac is not None:
        if samples_to_keep:
            max_mac = len(samples_to_keep) - max_mac
        elif samples_to_remove:
            max_mac = len(variations.samples) - len(samples_to_remove) - max_mac
        else:
            max_mac = len(variations.samples) - max_mac
        flt = MacFilter(max_mac=max_mac,
                        do_histogram=do_histogram,
                        do_filtering=True)
        pipeline.append(flt, id_='mac')

    if max_maf is not None:
        flt = MafFilter(max_maf=max_maf,
                        do_histogram=do_histogram,
                        do_filtering=True)
        pipeline.append(flt, id_='maf')

    if min_called is not None:
        flt = MinCalledGTsFilter(min_called=min_called,
                                 do_histogram=do_histogram,
                                 do_filtering=True)
        pipeline.append(flt, id_='called_rate')

    if max_het is not None:
        flt = ObsHetFilter(max_het=max_het,
                           min_call_dp=min_call_for_het,
                           do_histogram=do_histogram,
                           do_filtering=True)
        pipeline.append(flt, id_='obs_het')

    if out_variations:
        tier2_vars = out_variations
    elif cache_dir:
        tier2_vars = VariationsH5(str(h5_cache_path), 'w')
    else:
        tier2_vars = VariationsArrays()

    try:
        result = pipeline.run(variations, tier2_vars, chunk_size=chunk_size,
                            max_chunks_to_process=max_chunks_to_process,
                            kept_fields=kept_fields)
    except Exception:
        if cache_dir and h5_cache_path.exists():
            os.remove(h5_cache_path)
        raise

    if verbose:
        for step, res in result.items():
            if 'flt_stats' in res:
                stats = res['flt_stats']
                print(step)
                print('\tSNPs processed: %d' % stats['tot'])
                print('\tSNPs kept: %d' % stats['n_kept'])
                print('\tSNPs filtered_out: %d\n' % stats['n_filtered_out'])

    if do_histogram:
        hist_param_map = {'mac': {'hrange': None,
                                  'vlines': max_mac},
                          'called_rate': {'hrange': None,
                                          'vlines': min_called},
                          'obs_het': {'hrange': None,
                                      'vlines': max_het}}
        for step, res in result.items():
            if step in hist_param_map.keys():
                path = hist_path / (step + '.png')
                render_histogram(res, path.open('wb'), hist_param_map[step]['vlines'],
                                 hrange=hist_param_map[step]['hrange'])

    return tier2_vars


def get_low_quality_samples(variations, min_called_rate, sampling_rate=None,
                            regions_to_remove=None, n_bins=40,
                            missing_rate_hist_path=None, chunk_size=None,
                            max_chunks_to_process=None, out_dir=None):
    print('Getting low quality samples')
    if sampling_rate is not None or regions_to_remove is not None:
        kept_fields = [GT_FIELD]
        if regions_to_remove:
            kept_fields.extend([CHROM_FIELD, POS_FIELD])
        sampled_vars = VariationsArrays()
        filter_variations(variations, regions_to_remove=regions_to_remove,
                          sampling_rate=sampling_rate,
                          out_variations=sampled_vars,
                          chunk_size=chunk_size,
                          max_chunks_to_process=max_chunks_to_process,
                          kept_fields=kept_fields)
        print('Num variations to analyze:', sampled_vars.num_variations)
    else:
        sampled_vars = variations

    do_histogram = bool(missing_rate_hist_path)
    out_vars = VariationsArrays()
    res = filter_samples_by_missing_rate(sampled_vars,
                                         out_vars=out_vars,
                                         samples=None,
                                         min_called_rate=min_called_rate,
                                         n_bins=n_bins,
                                         do_histogram=do_histogram)

    if do_histogram:
        render_histogram(res, missing_rate_hist_path.open('wb'), min_called_rate)
    
    orig_samples = set(variations.samples)
    flt_samples = set(out_vars.samples)
    low_quality_samples = orig_samples.difference(flt_samples)
    return {'low_quality_samples': low_quality_samples}


def generate_blocks(variations, difference_rate_allowed=0.05,
                    min_num_gts_compared=10, chunk_size=100,
                    max_missing_rate_in_ref_snp=0.9, debug=False,
                    cache_dir=None):
    if cache_dir and not debug:
        key = ','.join(sorted(variations.samples)) + str(variations.num_variations)
        key += str(difference_rate_allowed)
        key += str(min_num_gts_compared)
        key += str(chunk_size)
        key += str(max_missing_rate_in_ref_snp)
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('blocks_' + key)
        if cache_path.exists():
            blocks = pickle.load(cache_path.open('rb'))
            return blocks
    else:
        cache_path = None

    print('Generating blocks')
    blocks = variation_generate_blocks(variations,
                                       difference_rate_allowed=difference_rate_allowed,
                                       min_num_gts_compared=min_num_gts_compared,
                                       chunk_size=chunk_size,
                                       max_missing_rate_in_ref_snp=max_missing_rate_in_ref_snp,
                                       debug=debug)
    if cache_dir:
        blocks = list(blocks)
        pickle.dump(blocks, cache_path.open('wb'))
    return blocks


def _select_var_with_lowest_missing_gts(variations):
    missing_rate_for_vars = calc_called_gt(variations)
    lowest_missing_rate_idx = numpy.argmax(missing_rate_for_vars)
    selected_vars = variations.get_chunk(slice(lowest_missing_rate_idx,
                                               lowest_missing_rate_idx + 1))
    return selected_vars


def keep_the_var_with_lowest_missing_gts_per_haplo_block(variations, cache_dir=None,
                                                         difference_rate_allowed=0.05):

    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += str(variations.num_variations)
        key += str(difference_rate_allowed)
        key = hashlib.md5(key.encode()).hexdigest()
        h5_cache_path = cache_dir / ('purged_vars' + key + '.h5')
        if h5_cache_path.exists():
            return VariationsH5(str(h5_cache_path), 'r')

    blocks = generate_blocks(variations,
                             difference_rate_allowed=difference_rate_allowed,
                             cache_dir=cache_dir)

    if cache_dir:
        out_variations = VariationsH5(str(h5_cache_path), 'w')
    else:
        out_variations = VariationsArrays()

    for block in blocks:
        variations_in_block = variations.get_chunk(slice(block['start_idx'],
                                                         block['stop_idx']))
        gts = _select_var_with_lowest_missing_gts(variations_in_block)[GT_FIELD]
        chunk = VariationsArrays()
        chunk[GT_FIELD] = gts
        chunk.samples = variations.samples[:]
        chunk[CHROM_FIELD] = numpy.array([block['chrom']])
        chunk[POS_FIELD] = numpy.array([block['start']])
        out_variations.put_chunks([chunk])

    return out_variations


def keep_most_variable_snp_per_window(variations, win_size, cache_dir=None):

    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += str(variations.num_variations)
        key += str(win_size)
        key = hashlib.md5(key.encode()).hexdigest()
        h5_cache_path = cache_dir / ('per_win_vars' + key + '.h5')
        if h5_cache_path.exists():
            return VariationsH5(str(h5_cache_path), 'r')

    if cache_dir:
        out_variations = VariationsH5(str(h5_cache_path), 'w')
    else:
        out_variations = VariationsArrays()

    for vars_in_win in variations.iterate_wins(win_size):
        if not vars_in_win.num_variations:
            continue
        elif vars_in_win == 1:
            chunk = vars_in_win
        else:
            mafs = calc_maf(vars_in_win, min_num_genotypes=0)
            max_max_snp_idx = numpy.argmin(mafs)
            chunk = vars_in_win.get_chunk(slice(max_max_snp_idx, max_max_snp_idx + 1))
        out_variations.put_chunks([chunk])

    return out_variations


def keep_less_missing_snp_per_window(variations, win_size, cache_dir=None):

    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += str(variations.num_variations)
        key += str(win_size)
        key = hashlib.md5(key.encode()).hexdigest()
        h5_cache_path = cache_dir / ('per_win_vars_less_missing' + key + '.h5')
        if h5_cache_path.exists():
            return VariationsH5(str(h5_cache_path), 'r')

    if cache_dir:
        out_variations = VariationsH5(str(h5_cache_path), 'w')
    else:
        out_variations = VariationsArrays()

    for vars_in_win in variations.iterate_wins(win_size):
        if not vars_in_win.num_variations:
            continue
        elif vars_in_win == 1:
            chunk = vars_in_win
        else:
            missing = calc_missing_gt(vars_in_win)
            max_snp_idx = numpy.argmin(missing)
            chunk = vars_in_win.get_chunk(slice(max_snp_idx, max_snp_idx + 1))
        out_variations.put_chunks([chunk])

    return out_variations


def keep_one_random_snp_per_window(variations, win_size, cache_dir=None):

    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += str(variations.num_variations)
        key += str(win_size)
        key = hashlib.md5(key.encode()).hexdigest()
        h5_cache_path = cache_dir / ('per_win_one_random' + key + '.h5')
        if h5_cache_path.exists():
            return VariationsH5(str(h5_cache_path), 'r')

    if cache_dir:
        out_variations = VariationsH5(str(h5_cache_path), 'w')
    else:
        out_variations = VariationsArrays()

    for vars_in_win in variations.iterate_wins(win_size):
        if not vars_in_win.num_variations:
            continue
        elif vars_in_win == 1:
            chunk = vars_in_win
        else:
            snp_idx = random.randrange(vars_in_win.num_variations)
            chunk = vars_in_win.get_chunk(slice(snp_idx, snp_idx + 1))
        out_variations.put_chunks([chunk])

    return out_variations
