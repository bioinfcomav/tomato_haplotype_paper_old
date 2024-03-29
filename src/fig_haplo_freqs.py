

# gene más introgresado en ecu
# gene más introgresado en pe
# algunos genes curiosos

import config

import pandas
import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from variation.variations import VariationsH5

import passport
import pop_building
import allele_freq
import colors
import matplotlib_support
import genome
import relevant_genes
import labels
from introgressions_for_genes import get_genes_with_most_introgressions


def get_gene_with_most_introgressions(founder_pop, target_pop, introgression_source_pop, pops):
    introgession_config = {'samples_in_pop_with_introgressions': sorted(pops[target_pop]),
                           'samples_in_founder_pop': sorted(pops[founder_pop]),
                           'samples_in_introgression_source_pop': sorted(pops[introgression_source_pop]),
                           'freq_threshold_to_consider_allele_present_in_founder_pop' : config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_PRESENT,
                           'freq_threshold_to_consider_allele_common_in_introgression_source_pop': config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_COMMON,
                            }

    introgression_freqs = get_genes_with_most_introgressions(variations,
                                                             introgession_config=introgession_config,
                                                             genes=genes, num_genes=100,
                                                             upstream_region=config.UPSTREAM_REGION_IN_BP,
                                                             cache_dir=config.CACHE_DIR)

    return list(introgression_freqs.keys())[0]


def _name_allele_freqs(freqs):
    most_freq_haplo_idx = numpy.argmax(freqs.loc[:, 'sll_mx'].values)
    new_index = list(freqs.index)
    new_index[most_freq_haplo_idx] = 'sll_mx'

    most_freq_haplo_idx2 = numpy.argmax(freqs.loc[:, 'sp_pe'].values)
    if most_freq_haplo_idx != most_freq_haplo_idx2:
        new_index[most_freq_haplo_idx2] = 'sp_pe'

    sp_ec_freqs = freqs.loc[:, 'sp_ec'].to_dict()
    allele_idxs_sorted_by_sp_ec_freq = reversed(sorted(enumerate(sp_ec_freqs.values()), key=lambda x: x[1]))

    for allele_idx, _ in allele_idxs_sorted_by_sp_ec_freq:
        if allele_idx not in [most_freq_haplo_idx, most_freq_haplo_idx2]:
            most_freq_haplo_idx3 = allele_idx
            break
    new_index[most_freq_haplo_idx3] = 'sp_ec'

    freqs.index = new_index

    return freqs


def _set_nice_labels(freqs):
    freqs.columns = [labels.LABELS[pop] for pop in  freqs.columns]
    return freqs

if __name__ == '__main__':

    debug = False

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, passports)

    ref_pops = ['sp_pe', 'sp_ec' ,'sll_mx']

    pop_order = ['sp_pe', 'sp_montane', 'sp_ec',
                 'slc_co',
                 'slc_ma', 'slc_ec', 'slc_pe', 'slc_world', 'sll_mx']

    if debug:
        n_succesful_attempts = 2
    else:
        n_succesful_attempts = 100
    res = allele_freq.calc_mean_haplo_allele_freqs(variations, ref_pops, pops,
                                                    n_succesful_attempts=n_succesful_attempts)
    typical_freqs = res['mean_freqs']

    pop_order = [labels.LABELS[pop] for pop in pop_order]

    genes = genome.Genes()

    founder_pop = 'slc_ma'

    target_pop = 'slc_pe'
    introgression_source_pop = 'sp_pe'
    most_introgressed_pe_gene_id = get_gene_with_most_introgressions(founder_pop, target_pop,
                                                                     introgression_source_pop, pops)

    target_pop = 'slc_ec'
    introgression_source_pop = 'sp_ec'
    most_introgressed_ec_gene_id = get_gene_with_most_introgressions(founder_pop, target_pop,
                                                                     introgression_source_pop, pops)

    haplo_colors = {'sp_pe': colors.HAPLO_COLORS['sp_peru'],
                    'sp_ec': colors.HAPLO_COLORS['sp_ecu'],
                    'sll_mx': colors.HAPLO_COLORS['sl'],
                    'other_1': colors.LIGHT_GRAY}
    color_schema = colors.ColorSchema(haplo_colors)

    fig = Figure((5, 5))
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes1 = fig.add_subplot(311)

    typical_freqs = _set_nice_labels(typical_freqs)
    allele_freq.plot_allele_freqs_bars(typical_freqs, axes1, pop_order=pop_order,
                                       color_schema=color_schema)
    chunk = allele_freq.get_chunk_for_gene(variations, genes,
                                           most_introgressed_ec_gene_id,
                                           upstream_region=config.UPSTREAM_REGION_IN_BP)
    freqs = allele_freq.calc_allele_freq(chunk, pops)
    freqs = _set_nice_labels(_name_allele_freqs(freqs))
    axes2 = fig.add_subplot(312)
    allele_freq.plot_allele_freqs_bars(freqs, axes2, pop_order=pop_order,
                                       color_schema=color_schema)


    chunk = allele_freq.get_chunk_for_gene(variations, genes, most_introgressed_pe_gene_id,
                                           upstream_region=config.UPSTREAM_REGION_IN_BP)

    freqs = allele_freq.calc_allele_freq(chunk, pops)
    freqs = _set_nice_labels(_name_allele_freqs(freqs))
    axes3 = fig.add_subplot(313)
    allele_freq.plot_allele_freqs_bars(freqs, axes3, pop_order=pop_order,
                                       color_schema=color_schema)

    if False:
        gene = {'common_name': 'edh1_mut',
                    'chrom': 'SL4.0ch09',
                    'start': 63111577,
                    'end': 63111579,
                    'min_num_var_for_haplo_freqs': 2}
        chunk = relevant_genes.get_variations_in_region(variations,
                                        chrom=gene['chrom'],
                                        start=gene['start'],
                                        end=gene['end'],
                                        min_num_vars=gene['min_num_var_for_haplo_freqs'])
        freqs = allele_freq.calc_allele_freq(chunk, pops)
        freqs = _set_nice_labels(_name_allele_freqs(freqs))
        axes4 = fig.add_subplot(414)
        allele_freq.plot_allele_freqs_bars(freqs, axes4, pop_order=pop_order,
                                        color_schema=color_schema)

    matplotlib_support.turn_off_x_axis(axes1)
    matplotlib_support.turn_off_x_axis(axes2)
    #matplotlib_support.turn_off_x_axis(axes3)

    matplotlib_support.set_axes_background(axes1)
    matplotlib_support.set_axes_background(axes2)
    matplotlib_support.set_axes_background(axes3)
    #matplotlib_support.set_axes_background(axes4)

    axes1.set_ylabel('Haplo. freq.')
    axes2.set_ylabel('Haplo. freq.')
    axes3.set_ylabel('Haplo. freq.')
    #axes4.set_ylabel('Haplo. freq.')

    axes1.set_title('Random region', fontsize=10)
    axes2.set_title(f'Gene with most introgressions in SLC Ec ({most_introgressed_ec_gene_id})', fontsize=10)
    axes3.set_title(f'Gene with most introgressions in SLC Pe ({most_introgressed_pe_gene_id})', fontsize=10)
    #axes4.set_title('Solyc09g075080 (EDH1)', fontsize=10)

    plot_path = config.FIG_HAPLO_FREQS
    fig.tight_layout()
    fig.savefig(plot_path)
