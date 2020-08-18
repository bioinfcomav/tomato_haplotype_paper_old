
import config

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5
from variation import CHROM_FIELD

import passport
from pop_building import get_pops
import matplotlib_support
from sample_selection import get_samples_for_criteria
from snp_filtering import filter_variations
from pcas_do import prepare_vars
from pop_distances import calc_kosman_dists
from pca import do_pcoa_from_dists
import colors
import labels
import fig_style


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
        axes.scatter(projections[mask,0], projections[mask,1], label=labels.LABELS[pop],
                     color=color, alpha=alpha, s=marker_size)

    axes.set_xlabel(f'Dim 1 ({var_percentages[0]:.1f}%)', fontsize=fig_style.X_LABEL_SMALL_FONT_SIZE)
    axes.set_ylabel(f'Dim 2 ({var_percentages[1]:.1f}%)', fontsize=fig_style.X_LABEL_SMALL_FONT_SIZE)
    return {'uniq_pops': uniq_pops}


def plot_samples_pca(axes, variations, criteria, passports, color_schema,
                     alpha, marker_size):
    vars_for_dists = get_vars_for_dists(variations, criteria, passports)

    dist_result = calc_kosman_dists(vars_for_dists, cache_dir=config.CACHE_DIR)
    multivar_result = do_pcoa_from_dists(dist_result)

    res = plot_pca(axes, multivar_result, pops_for_samples, color_schema, alpha=alpha, marker_size=marker_size)

    return res


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
    color_schema = colors.ColorSchema(colors.CLASSIFICATION_RANK1_COLORS)
    alpha = 0.7
    marker_size = 70


    if debug:
        max_chunks_to_process = 1000
        cache_dir = None
    else:
        max_chunks_to_process = None
        cache_dir = config.CACHE_DIR

    vars_path = config.WORKING_H5
    variations = VariationsH5(str(vars_path), 'r')

    sample_passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, sample_passports)
    pops_for_samples = {sample: pop for pop, samples in pops.items() for sample in samples}

    plot_path = config.FIG_ACC_PCA_HIERALCHICAL

    fig = Figure((10, 30))
    FigureCanvas(fig) # Don't remove it or savefig will fail later

    ['sp_montane', 'sll_vint', 'slc_world', None, 'slc_co', 'sll_mx', 'slc_ec', 'sp_ec', 'sll_modern',
    'sp_x_sp', 'sp_x_sl', 'slc_ma', 'slc_pe', 'sp_pe']

    popss = [['sp_pe', 'sp_montane', 'sp_ec', 'sp_x_sp'],
             ['slc_co', 'slc_pe_n', 'slc_pe_s', 'slc_ec', 'slc_ma'],
             ['slc_pe_n', 'slc_pe_s'],
             ['sll_mx', 'slc_world', 'slc_ma', 'slc_pe_n'],
             ['sll_modern', 'sll_vint','sll_mx']]

    axes_row_heights = [1 / len(popss)] * len(popss)

    uniq_pops = None
    for pops_idx, pops in enumerate(popss):
        axes = matplotlib_support.add_axes(fig, row_idx=pops_idx,
                                           axes_row_heights=axes_row_heights,
                                           bottom_margin=0.1)
        criteria = {'criteria': [((config.RANK1, pops, config.KEEP))], 'samples_to_remove': [],
                    'samples_to_keep': []}
        res = plot_samples_pca(axes, variations, criteria, passports=sample_passports,
                               color_schema=color_schema,
                               alpha=alpha, marker_size=marker_size)
        matplotlib_support.set_axes_background(axes, 'white')
        axes.legend(prop={'size': 14})
        if uniq_pops is None:
            uniq_pops = set(res['uniq_pops'])
        else:
            uniq_pops.update(res['uniq_pops'])

    fig.savefig(plot_path)
