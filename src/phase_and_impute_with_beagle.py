import config

import tempfile
import subprocess
from pathlib import Path
import gzip
from pprint import pprint
from functools import partial

import numpy

from variation.variations import VariationsH5
from variation.gt_parsers.vcf import VCFParser
from variation.variations.filters import SampleFilter, FLT_VARS
from variation import GT_FIELD, DP_FIELD, SNPS_PER_CHUNK, MISSING_INT, CHROM_FIELD
from variation.gt_writers.vcf import write_vcf
from variation.variations.pipeline import Pipeline

from solcap import get_solcap_markers, collect_locs_per_chrom

import check_imputation

BEAGLE_JAR = '/home/jope/soft/beagle.12Jul19.0df.jar'
BEAGLE_MEM = 100
DEFAULT_NE = 100000


def open_vcf(vcf_path, vcf_is_gzipped):
    if vcf_is_gzipped:
        fhand = gzip.open(vcf_path, 'rb')
    elif vcf_is_gzipped == False:
        fhand = vcf_path.open('rb')
    else:
        raise ValueError('Please indicate if vcf_is_gzipped')
    return fhand


def sort_vcf(in_unsorted_vcf_path, out_sorted_vcf_path):
    sort_process = subprocess.Popen(['vcf-sort', in_unsorted_vcf_path],
                                    stdout=subprocess.PIPE)
    stdout = out_sorted_vcf_path.open('wb')
    gzip_process = subprocess.Popen(['bgzip'],
                                    stdin=sort_process.stdout,
                                    stdout=stdout)
    sort_process.stdout.close()
    gzip_process.communicate()


def export_solcap_map(out_fhand):
    markers = get_solcap_markers(approx_phys_loc=True,
                                cache_dir=config.CACHE_DIR)
    genet_locs, phys_locs = collect_locs_per_chrom(markers)
    
    chroms = sorted(genet_locs.keys())
    for  chrom in chroms:
        for genet_pos, phys_pos in zip(genet_locs[chrom], phys_locs[chrom]):
            line = chrom + '\t.\t' + str(genet_pos) + '\t' + str(phys_pos) + '\n'
            out_fhand.write(line)


def _phase_and_impute_vcf_with_beagle(vcf_path, beagle_out_base_path,
                                      n_processes, stdout_fhand=None,
                                      stderr_fhand=None,
                                      samples_to_exclude=None,
                                      chroms_to_ignore=None,
                                      vcf_is_gzipped=None,
                                      ap=True, gp=True, reuse_files=False):

    cmd = ['java']
    cmd.append(f'-Xmx{BEAGLE_MEM}g')
    cmd.extend(['-jar', BEAGLE_JAR])
    cmd.append(f'map={config.BEAGLE_MAP}')
    cmd.append(f'ne={DEFAULT_NE}')
    cmd.append(f'gt={vcf_path}')
    cmd.append(f'out={beagle_out_base_path}')
    cmd.append(f'nthreads={n_processes}')
    cmd.append(f'ap={str(ap).lower()}')
    cmd.append(f'gp={str(gp).lower()}')

    if samples_to_exclude:
        fhand = tempfile.NamedTemporaryFile(suffix='.txt', mode='wt')
        fhand.write('\n'.join(samples_to_exclude))
        cmd.append(f'excludesamples={fhand.name}')
        fhand.flush()

    if chroms_to_ignore:
        chroms_fhand = tempfile.NamedTemporaryFile(prefix='beagle_exclude_markers.', mode='wt')
        chroms_to_ignore = [chrom.encode() for chrom in chroms_to_ignore]
        vcf_fhand = open_vcf(vcf_path, vcf_is_gzipped)
        for line in vcf_fhand:
            if line.startswith(b'#'):
                continue
            items = line.split()
            chrom, pos = items[0], items[1]
            if chrom not in chroms_to_ignore:
                continue
            chroms_fhand.write(f'{chrom}:{pos}\n')
        cmd.append(f'excludemarkers={chroms_fhand.name}')
        chroms_fhand.flush()

    phased_and_imputed_vcf = Path(str(beagle_out_base_path) + '.vcf.gz')
    if not reuse_files or not phased_and_imputed_vcf.exists():
        print('Running')
        print(' '.join(cmd))
        subprocess.run(cmd, stdout=stdout_fhand, stderr=stderr_fhand)

    if chroms_to_ignore:
        del chroms_fhand
    if samples_to_exclude:
        del fhand

    assert phased_and_imputed_vcf.exists()

    sorted_phased_and_imputed_vcf = Path(str(beagle_out_base_path) + '.sorted.vcf.gz')
    if not reuse_files or not sorted_phased_and_imputed_vcf.exists():
        sort_vcf(phased_and_imputed_vcf, sorted_phased_and_imputed_vcf)
    phased_and_imputed_vcf = sorted_phased_and_imputed_vcf
    phased_and_imputed_h5 = Path(str(beagle_out_base_path) + '.sorted.h5')

    if not reuse_files or not phased_and_imputed_h5.exists():
        print('Creating H5')
        fhand = gzip.open(phased_and_imputed_vcf, 'rb')
        vcf_parser = VCFParser(fhand=fhand)
        imputed_vars = VariationsH5(str(phased_and_imputed_h5), 'w')
        imputed_vars.put_vars(vcf_parser)


def phase_and_impute_vcf_with_beagle(vcf_path, beagle_out_base_path,
                                     n_processes, stdout_fhand=None,
                                     stderr_fhand=None,
                                     samples_to_exclude=None,
                                     chroms_to_ignore=None,
                                     vcf_is_gzipped=None,
                                     reuse_files=False):

    _phase_and_impute_vcf_with_beagle(vcf_path, beagle_out_base_path,
                                      n_processes,
                                      stdout_fhand=stdout_fhand,
                                      stderr_fhand=stderr_fhand,
                                      samples_to_exclude=samples_to_exclude,
                                      chroms_to_ignore=chroms_to_ignore,
                                      vcf_is_gzipped=vcf_is_gzipped,
                                      reuse_files=reuse_files)


def create_phased_h5(phased_vars, final_phased_and_imputed_vars,
                     phased_and_imputed_vars, orig_vars,
                     samples_to_keep_in_phased_and_imputed=None):

    phased_and_imputed_chunks = phased_and_imputed_vars.iterate_chunks(chunk_size=SNPS_PER_CHUNK)
    orig_chunks = orig_vars.iterate_chunks(chunk_size=SNPS_PER_CHUNK)

    if samples_to_keep_in_phased_and_imputed is not None:
        sample_filter = SampleFilter(samples=samples_to_keep_in_phased_and_imputed)

    for phased_and_imputed_chunk, orig_chunk in zip(phased_and_imputed_chunks, orig_chunks):
        phased_chunk = phased_and_imputed_chunk.copy()
        orig_chunk = orig_chunk.copy()
        if samples_to_keep_in_phased_and_imputed is not None:
            phased_chunk = sample_filter(phased_chunk)[FLT_VARS]
            orig_chunk = sample_filter(orig_chunk)[FLT_VARS]
        phased_chunk.metadata = orig_chunk.metadata
        gts = phased_chunk[GT_FIELD]
        orig_gts = orig_chunk[GT_FIELD]
        gts[orig_gts == MISSING_INT] = MISSING_INT
        phased_chunk[DP_FIELD] = orig_chunk[DP_FIELD]
        phased_vars.put_chunks([phased_chunk])

        if samples_to_keep_in_phased_and_imputed is not None:
            final_phased_and_imputed_chunk = sample_filter(phased_and_imputed_chunk)[FLT_VARS]
        else:
            final_phased_and_imputed_chunk = phased_and_imputed_chunk
        #print(final_phased_and_imputed_chunk.num_variations, len(final_phased_and_imputed_chunk.samples))
        final_phased_and_imputed_vars.put_chunks([final_phased_and_imputed_chunk])

    print('orig_vars:', orig_vars.num_variations, len(orig_vars.samples))
    print('final_phased_and_imputed_vars:', final_phased_and_imputed_vars.num_variations, len(final_phased_and_imputed_vars.samples))


def _accession_for_sample_is_desired(sample_id, passports, desired_accessions):
    passport = passports[sample_id]
    accession_id = passport.get('accession')
    if not accession_id:
        return False
    return accession_id in desired_accessions


if __name__ == '__main__':
    imputation_dir = config.IMPUTATION_DIR
    imputation_dir.mkdir(exist_ok=True)
    input_h5 = config.WORKING_H5
    phased_h5 = config.WORKING_PHASED_H5
    phased_and_imputed_h5 = config.WORKING_PHASED_AND_IMPUTED_H5
    vcf_path = config.WORKING_VCF

    create_vcf = False
    vcf_is_gzipped = False
    phase_and_impute = True
    check_imputation = False

    if create_vcf:
        print('Creating VCF')
        vars = VariationsH5(str(input_h5), 'r')
        fhand = gzip.open(vcf_path, 'wb')
        write_vcf(vars, out_fhand=fhand)
        fhand.close()
    if vcf_is_gzipped:
        beagle_out_base_path = imputation_dir / Path(vcf_path.stem).stem
    else:
        beagle_out_base_path = imputation_dir / Path(vcf_path.stem)

    stdout = imputation_dir / 'beagle.stdout'
    stderr = imputation_dir / 'beagle.stderr'
    if phase_and_impute:
        map_fhand = open(config.BEAGLE_MAP, 'w')
        export_solcap_map(out_fhand)
        phase_and_impute_vcf_with_beagle(vcf_path, beagle_out_base_path,
                                         stdout_fhand=stdout.open('wt'),
                                         stderr_fhand=stderr.open('wt'),
                                         chroms_to_ignore=['SL2.50ch00'],
                                         vcf_is_gzipped=vcf_is_gzipped,
                                         n_processes=1,
                                         reuse_files=False)

    if check_imputation:
        vars_for_check = VariationsH5(str(input_h5), 'r')
        check_result = check_imputation.check_imputation(vars_for_check, check_dir,
                                                         gt_rate_to_set_missing=0.1,
                                                         delete_tmp_files=False,
                                                         n_beagle_processes=10)

        # which samples can be imputed without too much error?
        samples = check_result['samples']
        fails = check_result['rate_failed_imputations_per_sample']
        samples_imputed_ok = [sample for sample, fail_rate in zip(samples, fails) if fail_rate < config.MAX_IMPUTATION_FAIL_RATE]
    else:
        beagle_phased_and_imputed_h5 = Path(str(beagle_out_base_path) + '.sorted.h5')
        phased_and_imputed_vars = VariationsH5(str(beagle_phased_and_imputed_h5), 'r')
        samples_imputed_ok = phased_and_imputed_vars.samples

    orig_vars = VariationsH5(str(input_h5), 'r')
    final_phased_and_imputed_vars = VariationsH5(str(phased_and_imputed_h5), 'w')
    phased_vars = VariationsH5(str(phased_h5), 'w')
    print(orig_vars.num_variations)
    print(phased_and_imputed_vars.num_variations)
    create_phased_h5(phased_vars, final_phased_and_imputed_vars,
                     phased_and_imputed_vars, orig_vars,
                     samples_imputed_ok)

    phased_vars = VariationsH5(str(phased_h5), 'r')
    chroms_orig = numpy.unique(orig_vars[CHROM_FIELD])
    chroms_phased = numpy.unique(phased_vars[CHROM_FIELD])
    assert not set(chroms_orig).difference(chroms_phased)
