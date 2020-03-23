
import os
from tempfile import NamedTemporaryFile
from pathlib import Path
import gzip
from pprint import pprint

import numpy

import config

from variation.variations.random import (copy_setting_gts_to_missing,
                                         sample_variations)
from variation.variations import VariationsH5
from variation.gt_writers.vcf import write_vcf
from variation.gt_parsers.vcf import VCFParser
from variation import SNPS_PER_CHUNK, GT_FIELD, MISSING_INT

from phase_and_impute_with_beagle import _phase_and_impute_vcf_with_beagle
from plot import plot_histogram, plot_scatter

DEBUG = False


def create_imputed_vars(in_vars, imputed_vars, vars_with_missing,
                        gt_rate_to_set_missing, var_sampling_rate,
                        delete_tmp_files, tmp_dir, n_beagle_processes):
    current_vars = vars
    orig_vars = vars
    # sampling
    if var_sampling_rate != 1:
        print('Sampling')
        sampled_fhand = NamedTemporaryFile(suffix='.sampled.h5',
                                           dir=tmp_dir,
                                           delete=delete_tmp_files)
        os.remove(sampled_fhand.name)
        out_vars = VariationsH5(sampled_fhand.name, mode='w')
        sample_variations(current_vars, out_vars=out_vars,
                          sample_rate=var_sampling_rate)
        current_vars = out_vars
        orig_vars = out_vars

    # set some gts to missing
    print('Setting GTs to missing')
    copy_setting_gts_to_missing(current_vars, out_vars=vars_with_missing,
                                gt_rate_to_missing=gt_rate_to_set_missing)

    #impute
    with_missing_base_fhand = NamedTemporaryFile(dir=tmp_dir,
                                                 suffix='.with_missing')
    with_missing_base_fpath = with_missing_base_fhand.name
    with_missing_base_fhand.close()
    with_missing_vcf_fpath = with_missing_base_fpath + '.vcf.gz'

    print('Writing VCF to impute')
    fhand = gzip.open(with_missing_vcf_fpath, 'wb')
    write_vcf(vars_with_missing, out_fhand=fhand)
    fhand.close()

    print('Imputing with Beagle')
    beagle_out_base_path = Path(with_missing_base_fpath + '.beagle.imputed')
    stdout_fhand = NamedTemporaryFile(suffix='.beagle.stdout', dir=tmp_dir,
                                      delete=delete_tmp_files)
    stderr_fhand = NamedTemporaryFile(suffix='.beagle.stderr', dir=tmp_dir,
                                      delete=delete_tmp_files)
    phase_and_impute_vcf_with_beagle(Path(with_missing_vcf_fpath),
                                     beagle_out_base_path,
                                     stdout_fhand=stdout_fhand,
                                     stderr_fhand=stderr_fhand,
                                     chroms_to_ignore=['SL2.50ch00'],
                                     vcf_is_gzipped=True,
                                     n_processes=n_beagle_processes)
    # generate_h5_from_imputed_vcf
    print('Calculating max field lens')
    imputed_vcf = beagle_out_base_path.with_suffix('.imputed.vcf.gz')

    fhand = gzip.open(imputed_vcf, 'rb')
    vcf_parser = VCFParser(fhand=fhand)
    for snp in vcf_parser.variations:
        continue
    max_field_lens = vcf_parser.max_field_lens
    max_field_str_lens = vcf_parser.max_field_str_lens

    print('Creating H5')
    fhand = gzip.open(imputed_vcf, 'rb')
    vcf_parser = VCFParser(fhand=fhand, max_field_lens=max_field_lens,
                           max_field_str_lens=max_field_str_lens)
    imputed_vars.put_vars(vcf_parser)

    print('Deleting tmp files')
    if var_sampling_rate != 1:
        sampled_fhand.close()
    if delete_tmp_files:
        if os.path.exists(with_missing_vcf_fpath):
            os.remove(with_missing_vcf_fpath)
        for path in beagle_out_base_path.glob('*'):
            os.remove(path)


def _find_out_if_gt_is_different(orig_gts, imputed_gts, gts_are_phased):
    if gts_are_phased:
        return numpy.all(orig_gts != imputed_gts, axis=2)
    else:
        are_equal1 = numpy.all(orig_gts == imputed_gts, axis=2)
        are_equal2 = numpy.logical_and(orig_gts[:, :, 0] == imputed_gts[:, :, 1],
                                       orig_gts[:, :, 1] == imputed_gts[:, :, 0])
        are_equal = numpy.logical_or(are_equal1, are_equal2)
        return numpy.logical_not(are_equal)


def compare_gts(orig_vars, with_more_missing_vars, imputed_vars):

    #  snp or sample missing rate vs errors
    orig_chunks = orig_vars.iterate_chunks(chunk_size=SNPS_PER_CHUNK)
    imputed_chunks = imputed_vars.iterate_chunks(chunk_size=SNPS_PER_CHUNK)
    with_more_missing_chunks = with_more_missing_vars.iterate_chunks(chunk_size=SNPS_PER_CHUNK)

    tot_wrong_imputations_per_sample = None
    tot_imputations_done_per_sample = None
    all_wrong_imputations_per_snp = None
    tot_missing_gts_per_sample = None
    for orig_chunk, with_more_missing_chunk, imputed_chunk in zip(orig_chunks, with_more_missing_chunks, imputed_chunks):
        orig_gts = orig_chunk[GT_FIELD]
        imputed_gts = imputed_chunk[GT_FIELD]
        with_more_missing_gts = with_more_missing_chunk[GT_FIELD]

        gt_was_set_to_missing = numpy.all(orig_gts != with_more_missing_gts,
                                          axis=2)
        #print(gt_was_set_to_missing)
        gt_is_different = _find_out_if_gt_is_different(orig_gts, imputed_gts,
                                                       gts_are_phased=False)
        #print(gt_is_different)
        gt_imputation_is_wrong = numpy.logical_and(gt_was_set_to_missing,
                                                   gt_is_different)
        wrong_imputations_per_sample = numpy.sum(gt_imputation_is_wrong, axis=0)
        wrong_imputations_per_snp = numpy.sum(gt_imputation_is_wrong, axis=1)
        gt_is_missing = numpy.any(orig_gts == MISSING_INT, axis=2)
        num_missing_gts_per_snp = numpy.sum(gt_is_missing, axis=1)
        num_missing_gts_per_sample = numpy.sum(gt_is_missing, axis=0)

        if tot_wrong_imputations_per_sample is None:
            tot_wrong_imputations_per_sample = wrong_imputations_per_sample
        else:
            tot_wrong_imputations_per_sample += wrong_imputations_per_sample

        if tot_imputations_done_per_sample is None:
            tot_imputations_done_per_sample = numpy.sum(gt_was_set_to_missing, axis=0)
        else:
            tot_imputations_done_per_sample += numpy.sum(gt_was_set_to_missing, axis=0)

        if all_wrong_imputations_per_snp is None:
            all_wrong_imputations_per_snp = wrong_imputations_per_snp
        else:
            all_wrong_imputations_per_snp = numpy.concatenate((all_wrong_imputations_per_snp,
                                                               wrong_imputations_per_snp))

        if tot_missing_gts_per_sample is None:
            tot_missing_gts_per_sample = num_missing_gts_per_sample
        else:
            tot_missing_gts_per_sample += num_missing_gts_per_sample

    #print(tot_wrong_imputations_per_sample)
    #print(tot_imputations_done_per_sample)
    #print(all_wrong_imputations_per_snp)
    #print(tot_missing_gts_per_sample)
    rate_failed_imputations_per_sample = tot_wrong_imputations_per_sample / tot_imputations_done_per_sample
    #print(rate_failed_imputations_per_sample * 100)
    return {'rate_failed_imputations_per_sample': rate_failed_imputations_per_sample,
            'rate_missing_gts_per_sample': tot_missing_gts_per_sample / orig_vars.num_variations,
            'samples': orig_vars.samples}


def check_imputation(vars, out_dir, gt_rate_to_set_missing,
                     var_sampling_rate=1, n_beagle_processes=1,
                     delete_tmp_files=True):
    out_dir.mkdir(exist_ok=True)

    tmp_dir = out_dir / 'tmp'
    tmp_dir.mkdir(exist_ok=True)

    if DEBUG:
        imputed_h5 = tmp_dir / 'tmpsu0t5e8f.imputed.h5'
        with_more_missing_h5 = tmp_dir / 'tmpbo6p0jii.with_some_missing_gts.h5'
    else:
        imputed_h5 = NamedTemporaryFile(suffix='.imputed.h5', dir=tmp_dir,
                                        delete=False)
        os.remove(imputed_h5.name)
        imputed_h5 = Path(imputed_h5.name)
        imputed_vars = VariationsH5(str(imputed_h5), 'w')
        with_more_missing_h5 = NamedTemporaryFile(suffix='.with_some_missing_gts.h5',
                                                  dir=tmp_dir, delete=False)
        os.remove(with_more_missing_h5.name)
        with_more_missing_h5 = Path(with_more_missing_h5.name)
        with_more_missing_vars = VariationsH5(str(with_more_missing_h5), 'w')
        create_imputed_vars(vars, imputed_vars, with_more_missing_vars,
                            gt_rate_to_set_missing, var_sampling_rate,
                            delete_tmp_files=False, tmp_dir=tmp_dir,
                            n_beagle_processes=n_beagle_processes)

    orig_vars = vars
    imputed_vars = VariationsH5(str(imputed_h5), 'r')
    vars_with_more_missing = VariationsH5(str(with_more_missing_h5), 'r')
    print(imputed_vars.num_variations)
    print(orig_vars.num_variations)
    print(vars_with_more_missing.num_variations)

    #compare per sample and per snp missing
    comparisons = compare_gts(orig_vars, vars_with_more_missing, imputed_vars)
    plot_path = out_dir / 'failed_imputations_per_sample_rate.hist.svg'
    plot_histogram(comparisons['rate_failed_imputations_per_sample'],
                   plot_path, labels={'x_axis': 'failed imputation rate',
                                      'y_axis': 'Num. samples'})

    plot_path = out_dir / 'failed_imputations_vs_missing_rate_per_sample_rate.hist.svg'
    plot_scatter(comparisons['rate_missing_gts_per_sample'],
                 comparisons['rate_failed_imputations_per_sample'],
                 plot_path, labels={'y_axis': 'failed imputation rate',
                                    'x_axis': 'Missing GT rate'})

    if delete_tmp_files:
        if imputed_h5.exists():
            os.remove(imputed_h5)
    return comparisons


if __name__ == '__main__':
    out_dir = config.IMPUTATION_CHECK_DIR
    vars = VariationsH5(str(config.TIER1_NO_LOW_QUAL_SAMPLES), mode='r')
    check_imputation(vars, out_dir, gt_rate_to_set_missing=0.1,
                     delete_tmp_files=False, n_beagle_processes=10)
