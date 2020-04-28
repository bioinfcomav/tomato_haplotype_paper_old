
import config

import os
import gzip
import subprocess

import numpy

import toyplot

from variation import GT_FIELD, MISSING_INT
from variation.variations import VariationsH5
from variation.matrix.stats import counts_by_row
from variation.variations.filters import SampleFilter, FLT_VARS

from passport import get_sample_passports
from pop_building import get_pops
from snp_filtering import (get_one_random_var_per_haplo_block,
                           keep_variations_variable_in_samples)
import trees
from plot import plot_scatter


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


def do_tree_mix_analysis_old(variations, out_dir, pops, num_migration_range, out_group=None,
                         n_boostraps=100):
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


def run_one_tree_mix_analysis(variations, out_dir, pops, num_migrations, out_group=None):

    snp_path = out_dir / 'snps.treemix.gz'
    snp_fhand = gzip.open(snp_path, 'wt')
    write_tree_mix_snp_file(variations, snp_fhand, pops)

    cmd = ['treemix', '-i', str(snp_path)]

    if out_group is not None:
        cmd.extend(['-root', out_group])

    out_base = out_dir / 'treemix_result'
    cmd.extend(['-o', str(out_base)])

    if num_migrations:
        cmd.extend(['-m', str(num_migrations)])

    stdout_fpath = str(out_base) + '.stdout'
    process = subprocess.run(cmd, check=True, stdout=open(stdout_fpath, 'wt'))

    newick_fpath = str(out_base) + '.treeout.gz'
    newick_str = gzip.open(newick_fpath, 'rt').read()

    tree_with_migrations = trees.parse_treemix_tree(newick_str)

    return {'tree': tree_with_migrations}


def do_tree_mix_analysis(variations, out_dir, pops, num_migration_range,
                         difference_rate_allowed_for_haplo_block,
                         max_maf=None,
                         blocks_cache_dir=None,
                         out_group=None, n_bootstraps=100):

    if n_bootstraps is None or n_bootstraps < 1:
        raise ValueError('n_bootstraps should be at least 1')

    for num_migrations in num_migration_range:
        this_iter_out_dir = out_dir / f'num_migrations_{num_migrations}'
        bootstrap_dir = this_iter_out_dir / 'boot'
        tree_list = trees.TreeWithMigrationsList()
        for boot_idx in range(n_bootstraps):
            this_iter_variations = get_one_random_var_per_haplo_block(variations,
                                                                      max_maf=max_maf,
                                                                      difference_rate_allowed=difference_rate_allowed_for_haplo_block,
                                                                      blocks_cache_dir=blocks_cache_dir)
            this_iter_iter_out_dir = bootstrap_dir / f'boot_{boot_idx}'
            os.makedirs(this_iter_iter_out_dir, exist_ok=True)
            res = run_one_tree_mix_analysis(this_iter_variations,
                                            out_dir=this_iter_iter_out_dir,
                                            pops=pops,
                                            num_migrations=num_migrations,
                                            out_group=out_group)
            tree_list.append(res['tree'])
        res = tree_list.get_consensus_tree()

        plot_path = this_iter_out_dir / f'majority_rule_consensus_tree.num_vars_{this_iter_variations.num_variations}.svg'
        draw_support = n_bootstraps > 1
        canvas = toyplot.Canvas(style={"background-color": "white"})
        axes = canvas.cartesian()
        trees.draw_tree_with_migrations(res['consensus_tree'], axes, draw_support=draw_support)
        axes.show = False
        toyplot.svg.render(canvas, str(plot_path))


def read_treemix_results(out_dir):
    num_migrationss = []
    likelihoods = []
    for path in out_dir.iterdir():
        if not path.is_dir():
            continue
        if not 'num_migrations' in path.name:
            continue

        num_migrations = int(path.name.split('_')[-1])

        base_boot_dir = path / 'boot'
        this_likelihoods = []
        for boot_dir in base_boot_dir.iterdir():
            if 'boot' not in boot_dir.name:
                continue
            treemix_stdout_path = boot_dir / 'treemix_result.stdout'
            if not treemix_stdout_path.exists():
                continue
            lines = treemix_stdout_path.open('rt').read().splitlines()
            lines = list(filter(lambda x: 'ln(likelihood)' in x, lines))
            likelihood = float(lines[-1].split(':')[-1].strip())
            this_likelihoods.append(likelihood)
        likelihood = numpy.mean(this_likelihoods)
        num_migrationss.append(num_migrations)
        likelihoods.append(likelihood)
    return {'num_migrations': num_migrationss, 'likelihoods': likelihoods}


def plot_likelihoods(out_dir):
    res = read_treemix_results(out_dir)

    plot_path = out_dir / 'likelihoods.svg'
    plot_scatter(res['num_migrations'], res['likelihoods'], plot_path,
                 labels={'x_axis': 'Num. migrations', 'y_axis': 'Likelihood'})



if __name__ == '__main__':

    difference_rate_allowed_for_haplo_block = 0.40
    max_maf = 0.95
    cache_dir = config.CACHE_DIR

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
        num_migration_range = range(0, 5)

    if True:
        main_pops = ['sp_pe' ,'sp_ec', 'slc_ma', 'slc_ec', 'slc_pe', 'sll_mx']
        all_pops = main_pops
        out_fname = 'only_america'
        out_group = 'sp_pe'
        num_migration_range = range(0, 8)

    pops_descriptions = {config.RANK1: all_pops}
    pops = get_pops(pops_descriptions, passports)
    all_samples = {sample for samples in pops.values() for sample in samples}

    variations = keep_variations_variable_in_samples(variations, all_samples)
    variations = SampleFilter(all_samples)(variations)[FLT_VARS]

    out_fname2 = f'max_ld_corr_{difference_rate_allowed_for_haplo_block}'
    if max_maf is not None:
        out_fname2 += f'.maf_maf_{max_maf}'
    out_dir = config.TREE_MIX_DIR / out_fname / out_fname2
    os.makedirs(out_dir, exist_ok=True)

    do_tree_mix_analysis(variations, out_dir, pops, num_migration_range=num_migration_range,
                        difference_rate_allowed_for_haplo_block=difference_rate_allowed_for_haplo_block,
                        out_group=out_group, n_bootstraps=n_boostraps,
                        max_maf=max_maf,
                        blocks_cache_dir=cache_dir)

    plot_likelihoods(out_dir)