
import config

from collections import defaultdict
from pprint import pprint
import time
import hashlib
import pickle
import gzip
import math

import numpy
from pandas import DataFrame

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib
from matplotlib.lines import Line2D

import cartopy.crs as ccrs

from variation import CHROM_FIELD
from variation.variations import VariationsH5
from variation.variations.filters import VarsSamplingFilter2, FLT_VARS
from variation.clustering import do_tree

from passport import get_sample_passports
from sample_selection import get_samples_for_criteria
from snp_filtering import (filter_variations,
                           keep_most_variable_snp_per_window,
                           keep_the_var_with_lowest_missing_gts_per_haplo_block)
import colors
from plot import plot_var_density_along_genome
from pop_distances import calc_kosman_dists
from pca import (do_pcoa_from_dists, write_multivariant_result_for_curly,
                 write_pca_curlywhirly_file)
from pcas_do import get_sample_selection_criteria, prepare_vars
import passport
from pop_building import get_pops
import trees
import labels
from plot_geo_map import plot_geo_rank1_for_main_pops
from haplo_auto_classification import detected_outliers_and_classify_haplos
from haplo_pca_plotting import filter_aligned_pcoas_df_for_samples
from fig1b import (LEGEND_FONT_SIZE, HAPLO_PCOAS_X_LIMS, HAPLO_PCOAS_Y_LIMS,
                   X_LABEL_SMALL_FONT_SIZE, LEGEND_MARKER_SIZE)


def get_vars_for_dists(variations, criteria, passports):
    all_samples = variations.samples

    samples_to_use = get_samples_for_criteria(all_samples,
                                              passports,
                                              criteria,
                                              skip_samples_with_no_passport=config.SKIP_SAMPLES_WITH_NO_PASSPORT)

    if filter_by_maf:
        max_maf = config.TIER2_PCA_MAX_MAF
    else:
        max_maf = None

    print(variations.num_variations)
    tier2_vars = filter_variations(variations,
                                   chunk_size=chunk_size,
                                   samples_to_keep=samples_to_use,
                                   cache_dir=cache_dir,
                                   filter_out_vars_with_non_major_allele_count_le=config.TIER2_MAX_MAC,
                                   max_maf=max_maf,
                                   min_called=config.TIER2_MIN_CALLED,
                                   max_het=config.TIER2_MAX_HET,
                                   min_call_for_het=config.TIER2_MAX_HET_MIN_CALL_DP,
                                   kept_fields=config.RELEVANT_FIELDS,
                                   max_chunks_to_process=max_chunks_to_process,
                                   remove_non_variable_snvs=True,
                                   verbose=True
                                   )

    if not tier2_vars.num_variations:
        raise ValueError('No SNPs left after tier2')

    print(f'tier2 vars: {tier2_vars.num_variations}')

    var_counts_per_chrom = dict(zip(*numpy.unique(tier2_vars[CHROM_FIELD], return_counts=True)))

    if not debug and len(var_counts_per_chrom) < 12:
        print('Only some chromosomes have SNPs:')
        for chrom, num_vars in var_counts_per_chrom.items():
            chrom = chrom.decode()
            print(f'\t{chrom}:\t{num_vars}')
        raise RuntimeError('Some chromosomes have no SNPs')

    vars_for_dists = prepare_vars(tier2_vars,
                                  plot_var_density=plot_var_density,
                                  win_size_to_keep_only_the_most_variable_var=win_size_to_keep_only_the_most_variable_var,
                                  keep_at_most_n_vars=keep_at_most_n_vars,
                                  difference_rate_allowed_for_haplo_block=difference_rate_allowed_for_haplo_block,
                                  min_num_snps_to_use=min_num_snps_to_use,
                                  ignore_chromosome_representation_check=ignore_chromosome_representation_check,
                                  cache_dir=cache_dir)
    return vars_for_dists


def plot_pca(axes, multivar_result, pops_for_samples, color_schema, alpha, marker_size):

    projections = multivar_result['projections']
    var_percentages = multivar_result['var_percentages']
    samples = multivar_result['samples']

    pops = numpy.array([pops_for_samples[sample] for sample in samples])
    uniq_pops = list(set(pops))

    for pop in uniq_pops:
        mask = pops == pop
        color = color_schema[pop]
        axes.scatter(projections[mask,0], projections[mask,1], label=labels.LABELS[pop], color=color, alpha=alpha, s=marker_size)

    axes.set_xlabel(f'Dim 1 ({var_percentages[0]:.1f}%)', fontsize=X_LABEL_SMALL_FONT_SIZE)
    axes.set_ylabel(f'Dim 2 ({var_percentages[1]:.1f}%)', fontsize=X_LABEL_SMALL_FONT_SIZE)


def plot_haplo_pca_for_samples(aligned_pcoas_df, haplo_classification,
                               samples, axes,
                               x_lims, y_lims,
                               color_schema=None,
                               alpha=1):
    aligned_pcoas_df_for_sample = filter_aligned_pcoas_df_for_samples(aligned_pcoas_df,
                                                                      samples)
    
    haplo_classes = numpy.array([haplo_classification[haplo_id] for haplo_id in aligned_pcoas_df_for_sample.index])
    uniq_haplo_classes = sorted(set(haplo_classes))
    
    for haplo_class in uniq_haplo_classes:
        mask = haplo_classes == haplo_class
        aligned_pcoas_for_sample_and_haplo_class = aligned_pcoas_df_for_sample.values[mask, :]
        color = colors.HAPLO_COLORS.get(haplo_class, None)
        axes.scatter(aligned_pcoas_for_sample_and_haplo_class[:, 0],
                     aligned_pcoas_for_sample_and_haplo_class[:, 1],
                     color=color, label=labels.HAPLO_LABELS[haplo_class],
                     alpha=alpha)


def create_figure(variations, pops, aligned_pcoas_df, haplo_classification, sample_for_haplo_pca):

    pops_for_samples = {sample: pop for pop, samples in pops.items() for sample in samples}

    alpha = 0.8
    marker_size = 40

    out_dir = config.PAPER_FIGURES_DIR
    out_dir.mkdir(exist_ok=True)
    plot_path = out_dir / 'fig0.png'

    fig = Figure((20, 20))
    FigureCanvas(fig) # Don't remove it or savefig will fail later

    ncols = 6
    nrows = 4
    height_ratios = [0.1, 1, 0.55, 0.55]
    wspace = 0
    hspace = 0
    grid_spec = matplotlib.gridspec.GridSpec(nrows, ncols, figure=fig,
                                             height_ratios=height_ratios,
                                             wspace=wspace, hspace=hspace)

    legend_axes = fig.add_subplot(grid_spec[0, :])
    plot_pop_legend(pops, legend_axes)

    # PCA with all
    criteria = {'criteria': [], 'samples_to_remove': [],
                'samples_to_keep': []}

    vars_for_dists = get_vars_for_dists(variations, criteria, passports)

    dist_result = calc_kosman_dists(vars_for_dists, cache_dir=config.CACHE_DIR)
    multivar_result = do_pcoa_from_dists(dist_result)

    color_schema = colors.ColorSchema(colors.CLASSIFICATION_RANK1_COLORS)
    axes1 = add_axes(fig)

    plot_pca(axes1, multivar_result, pops_for_samples, color_schema, alpha=alpha, marker_size=marker_size)
    fig.savefig(plot_path)
    return

    axes2 = fig.add_subplot(grid_spec[1, 2:4])
    pops_for_pca = ['sll_modern', 'sll_vint','sll_mx', 'sll_old_cultivars', 'slc_world', 'slc_ma']
    criteria = {'criteria': [((config.RANK1, pops_for_pca, config.KEEP))], 'samples_to_remove': [],
                'samples_to_keep': []}
    vars_for_dists = get_vars_for_dists(variations, criteria, passports)
    dist_result = calc_kosman_dists(vars_for_dists, cache_dir=config.CACHE_DIR)
    multivar_result = do_pcoa_from_dists(dist_result)
    plot_pca(axes2, multivar_result, pops_for_samples, color_schema, alpha=alpha, marker_size=marker_size)

    axes1.set_facecolor('white')
    axes2.set_facecolor('white')

    map_axes = fig.add_subplot(grid_spec[1, 4:], projection=ccrs.PlateCarree())
    plot_geo_rank1_for_main_pops(sample_passports, axes=map_axes, plot_legend=False)

    pop_order = ['sp_pe', 'sp_ec', 'slc_co', 'slc_ma', 'slc_world', 'sll_mx', 
                 'sp_pe_inter-andean', 'slc_ec', 'slc_pe', 'sll_vint', 'sll_modern', 'sp_x_sl']
                    #'sll_old_cultivars', , , None]
    for pop_idx, pop in enumerate(pop_order):
        col = pop_idx % ncols
        row = pop_idx  // ncols + 2
        pop_axes = fig.add_subplot(grid_spec[row, col])
        samples = pops[pop]
        plot_haplo_pca_for_samples(aligned_pcoas_df, haplo_classification, samples,
                                   pop_axes, HAPLO_PCOAS_X_LIMS, HAPLO_PCOAS_Y_LIMS,
                                   color_schema=None, alpha=0.05)
        pop_axes.set_facecolor('white')
        turn_off_both_axis(axes=pop_axes)
        pop_axes.set_xlim(HAPLO_PCOAS_X_LIMS)
        pop_axes.set_ylim(HAPLO_PCOAS_Y_LIMS)
        set_axis_color(pop_axes)
        
        pop_axes.text(-0.03, 0.03, labels.LABELS[pop], fontsize=14, horizontalalignment='center')

    fig.savefig(plot_path)


if __name__ == '__main__':

    debug = False

    filter_by_maf = True
    chunk_size = 1000
    plot_var_density = False
    win_size_to_keep_only_the_most_variable_var = 1e5
    keep_at_most_n_vars = 2000
    difference_rate_allowed_for_haplo_block = 0.2
    max_num_vars_for_pca = 2000
    min_num_snps_to_use = 100
    ignore_chromosome_representation_check = False
    sample_for_haplo_pca = 'cervil'

    if debug:
        max_chunks_to_process = 1000
        cache_dir = None
    else:
        max_chunks_to_process = None
        cache_dir=config.CACHE_DIR

    passports = get_sample_passports()

    vars_path = config.WORKING_H5
    variations = VariationsH5(str(vars_path), 'r')

    sample_passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, sample_passports)

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    imputed_variations = VariationsH5(str(vars_path), 'r')
    samples_to_use = {sample for samples in pops.values() for sample in samples}
    res = detected_outliers_and_classify_haplos(variations,
                                                win_params={'min_num_snp_for_window': config.MIN_NUM_SNPS_FOR_HAPLO_IN_PCA,
                                                            'win_size': config.HAPLO_WIN_SIZE},
                                                num_wins_to_process=None,
                                                samples_to_use=samples_to_use,
                                                n_dims_to_keep=config.N_DIMS_TO_KEEP,
                                                classification_config=config.CLASSIFICATION_CONFIG,
                                                classification_outlier_config=config.CLASSIFICATION_OUTLIER_CONFIG,
                                                outlier_configs=config.OUTLIER_CONFIGS,
                                                out_dir=config.HAPLO_PCOA_DIR,
                                                pops=pops,
                                                outliers_return_aligned_pcoas=False,
                                                only_outliers=False,
                                                classification_references=config.CLASSIFICATION_REFERENCES,
                                                supervised_classification_config=config.SUPERVISED_CLASSIFICATION_CONFIG,
                                                cache_dir=cache_dir)

    create_figure(variations, pops,
                  res['aligned_pcoas_df'], res['classification'], sample_for_haplo_pca=sample_for_haplo_pca)
