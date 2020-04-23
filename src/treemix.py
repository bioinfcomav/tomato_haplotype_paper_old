
import config

import os
import gzip
import subprocess

import numpy

from variation import GT_FIELD, MISSING_INT
from variation.variations import VariationsH5
from variation.matrix.stats import counts_by_row
from variation.variations.filters import SampleFilter, FLT_VARS

from passport import get_sample_passports
from pop_building import get_pops
from snp_filtering import (keep_the_var_with_lowest_missing_gts_per_haplo_block,
                           keep_variations_variable_in_samples)


def count_minor_and_major_allele_counts_per_pop(variations, pops):
    alleles = numpy.unique(variations[GT_FIELD])
    alleles = set(alleles).difference([MISSING_INT])
    alleles = sorted(alleles)

    num_vars = variations.num_variations

    counts = counts_by_row(variations[GT_FIELD], missing_value=MISSING_INT,
                            alleles=alleles)
    major_alleles_idx = numpy.argmax(counts, axis=1)

    all_samples = variations.samples
    for pop, samples in pops.items():
        variations_for_pop = SampleFilter(samples)(variations)[FLT_VARS]
        counts = counts_by_row(variations_for_pop[GT_FIELD],
                               missing_value=MISSING_INT,
                               alleles=alleles)
        vars_idx = numpy.arange(num_vars)
        major_allele_count = counts[vars_idx, major_alleles_idx]
        minor_allele_count = numpy.sum(counts, axis=1) - major_allele_count
        yield {'pop': pop,
               'major_allele_count': major_allele_count,
               'minor_allele_count': minor_allele_count}


def write_tree_mix_snp_file(variations, fhand, pops):

    num_vars = variations.num_variations

    alelle_counts = count_minor_and_major_allele_counts_per_pop(variations, pops)

    comma_sep = numpy.array([','] * num_vars).reshape(1, num_vars)
    space_sep = numpy.array([' '] * num_vars).reshape(1, num_vars)
    new_lines = numpy.array(['\n'] * num_vars).reshape(1, num_vars)

    pops = []
    to_write = None
    for counts in alelle_counts:
        pops.append(counts['pop'])

        if to_write is not None:
            to_write = numpy.char.add(to_write, space_sep)

        major_counts = numpy.array(counts['major_allele_count'], str).reshape(1, num_vars)
        if to_write is None:
            to_write = major_counts
        else:
            to_write = numpy.char.add(to_write, major_counts)

        to_write = numpy.char.add(to_write, comma_sep)

        to_write = numpy.char.add(to_write, numpy.array(counts['minor_allele_count'], str).reshape(1, num_vars))
    to_write = numpy.char.add(to_write, new_lines)

    to_write = to_write.squeeze()

    fhand.write(' '.join(pops))
    fhand.write('\n')

    for idx in range(to_write.shape[0]):
        fhand.write(to_write[idx])

    fhand.flush()


def do_tree_mix_analysis(variations, out_dir, pops, num_migration_range, out_group=None):
    tmp_dir = out_dir / 'tmp'
    tmp_dir.mkdir(exist_ok=True)

    snp_path = tmp_dir / 'snps.treemix.gz'
    snp_fhand = gzip.open(snp_path, 'wt')
    write_tree_mix_snp_file(variations, snp_fhand, pops)

    cmd = ['treemix', '-i', str(snp_path)]

    if out_group is not None:
        cmd.extend(['-root', out_group])

    for num_migrations in num_migration_range:
        this_iter_cmd = cmd[:]
    
        this_iter_out_dir = out_dir / f'num_migrations_{num_migrations}'
        this_iter_out_dir.mkdir(exist_ok=True)
        this_iter_out_base = this_iter_out_dir / 'treemix_result'
        this_iter_cmd.extend(['-o', str(this_iter_out_base)])

        if num_migrations:
            this_iter_cmd.extend(['-m', str(num_migrations)])
        #print(' '.join())
        subprocess.run(this_iter_cmd, check=True)

        tree_plot_fpath = str(this_iter_out_base) + '.tree.svg'
        r_plot_command = f'source("{config.TREE_MIX_PLOTTING_R_SOURCE}");pdf(NULL);svg("{tree_plot_fpath}");plot_tree("{this_iter_out_base}");dev.off()'
        subprocess.run(['R', '-e', r_plot_command])


if __name__ == '__main__':

    cache_dir = config.CACHE_DIR
    difference_rate_allowed_for_haplo_block = 0.5

    vars_path = config.WORKING_PHASED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = get_sample_passports()

    n_boostraps = 100

    if True:
        main_pops = ['sp_pe' ,'sp_ec', 'slc_ma', 'slc_ec', 'slc_pe', 'sll_mx']
        vintage_pops = ['sll_old_cultivars', 'sll_vint', 'slc_world']
        #hybrid_pops = ['sll_modern', 'sp_x_sl']
        all_pops = main_pops + vintage_pops #+ hybrid_pops
        out_fname = 'up_to_vintage'
        out_group = 'sp_pe'
        num_migration_range = range(0, 6)

    if True:
        main_pops = ['sp_pe' ,'sp_ec', 'slc_ma', 'slc_ec', 'slc_pe', 'sll_mx']
        all_pops = main_pops
        out_fname = 'only_america'
        out_group = 'sp_pe'
        num_migration_range = range(0, 6)

    pops_descriptions = {config.RANK1: all_pops}
    pops = get_pops(pops_descriptions, passports)
    all_samples = {sample for samples in pops.values() for sample in samples}

    variations = keep_variations_variable_in_samples(variations, all_samples)
    variations = SampleFilter(all_samples)(variations)[FLT_VARS]
    variations = keep_the_var_with_lowest_missing_gts_per_haplo_block(variations,
                                                                      difference_rate_allowed=difference_rate_allowed_for_haplo_block,
                                                                      cache_dir=cache_dir)

    out_dir = config.TREE_MIX_DIR / out_fname / f'max_ld_corr_{difference_rate_allowed_for_haplo_block}'
    os.makedirs(out_dir, exist_ok=True)
    
    do_tree_mix_analysis(variations, out_dir, pops, num_migration_range=num_migration_range,
                         out_group=out_group)
