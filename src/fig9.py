

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
import introgressions_for_genes
import introgressions
import genome

print('TODO fig 10:, genes relacionados con la domesticación')


def get_gene_with_most_introgressions(founder_pop, target_pop, introgression_source_pop, pops):
    introgession_config = {'samples_in_pop_with_introgressions': pops[target_pop],
                            'samples_in_founder_pop': pops[founder_pop],
                            'samples_in_introgression_source_pop': pops[introgression_source_pop],
                            'freq_threshold_to_consider_allele_present_in_founder_pop' : config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_PRESENT,
                            'freq_threshold_to_consider_allele_common_in_introgression_source_pop': config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_COMMON,
                            }
    introgression_freqs = introgressions.calc_introgession_freq_for_vars(variations, introgession_config)

    introgression_freqs = introgressions_for_genes.get_genes_with_most_introgressions(variations, introgression_freqs, genes, num_genes=1)
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

if __name__ == '__main__':

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, passports)

    ref_pops = ['sp_pe', 'sp_ec' ,'sll_mx']

    pop_order = ['sp_pe', 'sp_pe_inter-andean', 'sp_ec',
                 'slc_co',
                 'slc_ma', 'slc_ec', 'slc_pe', 'slc_world', 'sll_mx']

    if False:
        res = allele_freq.calc_mean_haplo_allele_freqs(variations, ref_pops, pops, n_succesful_attempts=100)
        typical_freqs = res['mean_freqs']
    else:
        pop_names = [None, 'slc_co', 'slc_ec', 'slc_ma', 'slc_pe', 'slc_world', 'sll_modern', 'sll_mx', 'sll_old_cultivars', 'sll_vint', 'sp_ec', 'sp_pe', 'sp_pe_inter-andean', 'sp_x_sl', 'sp_x_sp']
        freqs = numpy.array([[0.00041667, 0.        , 0.00297619, 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.0040625 , 0.00478261, 0.        , 0.00123077, 0.002     ],
       [0.00010417, 0.        , 0.0002381 , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.0059375 , 0.00434783, 0.01      , 0.00053846, 0.0035    ],
       [0.00947917, 0.01      , 0.00678571, 0.01      , 0.00987805,
        0.01      , 0.01      , 0.01      , 0.01      , 0.01      ,
        0.        , 0.00043478, 0.        , 0.00823077, 0.0045    ],
       [0.        , 0.        , 0.        , 0.        , 0.00012195,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.00043478, 0.        , 0.        , 0.        ]])
        typical_freqs= pandas.DataFrame(freqs,
                                        index=['sp_pe', 'sp_ec', 'sll_mx', 'other_1'],
                                        columns=pop_names)

    genes = genome.Genes()

    founder_pop = 'slc_ma'
    if False:
        target_pop = 'slc_pe'
        introgression_source_pop = 'sp_pe'
        most_introgressed_pe_gene_id = get_gene_with_most_introgressions(founder_pop, target_pop, introgression_source_pop, pops)
    else:
        most_introgressed_pe_gene_id = 'Solyc07g061780'

    if False:
        target_pop = 'slc_ec'
        introgression_source_pop = 'sp_ec'
        most_introgressed_ec_gene_id = get_gene_with_most_introgressions(founder_pop, target_pop, introgression_source_pop, pops)
    else:
        most_introgressed_ec_gene_id = 'Solyc06g062400'

    haplo_colors = {'sp_pe': colors.HAPLO_COLORS['sp_peru'],
                    'sp_ec': colors.HAPLO_COLORS['sp_ecu'],
                    'sll_mx': colors.HAPLO_COLORS['sl'],
                    'other_1': colors.LIGHT_GRAY}
    color_schema = colors.ColorSchema(haplo_colors)

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes1 = fig.add_subplot(311)

    allele_freq.plot_allele_freqs_bars(typical_freqs, axes1, pop_order=pop_order,
                                       color_schema=color_schema)
    matplotlib_support.set_axes_background(axes1)

    chunk = allele_freq.get_chunk_for_gene(variations, genes, most_introgressed_ec_gene_id)
    freqs = allele_freq.calc_allele_freq(chunk, pops)
    freqs = _name_allele_freqs(freqs)
    axes2 = fig.add_subplot(312)
    allele_freq.plot_allele_freqs_bars(freqs, axes2, pop_order=pop_order,
                                       color_schema=color_schema)


    chunk = allele_freq.get_chunk_for_gene(variations, genes, most_introgressed_pe_gene_id)
    freqs = allele_freq.calc_allele_freq(chunk, pops)
    freqs = _name_allele_freqs(freqs)
    axes3 = fig.add_subplot(313)
    allele_freq.plot_allele_freqs_bars(freqs, axes3, pop_order=pop_order,
                                       color_schema=color_schema)

    matplotlib_support.turn_off_x_axis(axes1)
    matplotlib_support.turn_off_x_axis(axes2)

    print('TODO edh1 with imputed vars')
    print('TODO'*100)
    print('Rerun everything with imputed vars')

    plot_path = config.FIGURES_DIR / 'fig9.svg'
    fig.tight_layout()
    fig.savefig(plot_path)

