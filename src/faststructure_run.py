
import config

import os
import hashlib
import gzip
import subprocess
from pathlib import Path
from random import randint

from variation.variations import VariationsH5
from variation.gt_writers.vcf import write_vcf

from snp_filtering import filter_variations


def _create_plink_variation_files(variations, plink_files_dir, if_exists_skip_creation=False):
    key = ','.join(sorted(variations.samples))
    key += str(variations.num_variations)
    key = hashlib.md5(key.encode()).hexdigest()

    base_path = plink_files_dir / (str(key) + '.plink')
    bim_fpath = str(base_path) + '.bim'
    if os.path.exists(bim_fpath) and if_exists_skip_creation:
        return {'base_path': base_path}

    vcf_path = plink_files_dir / (str(key) + '.vcf.gz')
    vcf_fhand = gzip.open(vcf_path, 'wb')
    write_vcf(variations, vcf_fhand)

    cmd = [config.PLINK_BIN, '--vcf', str(vcf_path),
           '--out',  str(base_path),
           '--allow-extra-chr', '--double-id', '--vcf-half-call', 'missing',
           '--set-missing-var-ids', '@:#', '--make-bed']
    subprocess.run(cmd, check=True)

    return {'base_path': base_path}


def _run_faststructure(k, input_fpath, output_fpath, prior):

    cmd = [str(config.PYTHON2_BIN_FOR_FASTSTRUCTURE),
           str(config.FASTSTRUCTURE_BIN)]
    cmd.extend(['-K', str(k)])
    cmd.extend(['--input', str(input_fpath)])
    cmd.extend(['--output', str(output_fpath)])
    cmd.extend(['--prior', prior])
    cmd.extend(['--seed', str(randint(100, 1000))])
    cmd.extend(['--full'])
    stderr_path = Path(str(output_fpath) + '.structure.{}.stderr'.format(k))
    stdout_path = Path(str(output_fpath) + '.structure.{}.stdout'.format(k))

    print('Running: ', ' '.join(cmd))
    subprocess.run(cmd, stdout=stdout_path.open('wt'),
                   stderr=stderr_path.open('wt'), check=True)


def run_faststructure(variations, k_range, prior='logistic'):
    out_dir = config.FASTSTRUCTURE_DIR / prior
    run_dir = config.FASTSTRUCTURE_RUN_DIR / prior

    plink_files_dir = config.FASTSTRUCTURE_PLINK_DIR
    os.makedirs(plink_files_dir, exist_ok=True)

    res = _create_plink_variation_files(variations, plink_files_dir, if_exists_skip_creation=True)
    plink_base_path = res['base_path']

    for k in range(*k_range):
        run_dir_for_this_k = run_dir / f'k_{k}'
        out_base_path = run_dir_for_this_k / config.FASTSTRUCTURE_RESULT_BASE_FNAME
        os.makedirs(run_dir_for_this_k, exist_ok=True)
        _run_faststructure(k, str(plink_base_path), str(out_base_path), prior)


if __name__ == '__main__':

    k_range = (2, 11)
    max_maf = 0.95

    vars_path = config.WORKING_H5
    variations = VariationsH5(str(vars_path), 'r')

    variations = filter_variations(variations, max_maf=max_maf)

    run_faststructure(variations, k_range, prior='simple')
