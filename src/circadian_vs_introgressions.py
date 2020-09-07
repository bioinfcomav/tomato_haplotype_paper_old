
import config

import subprocess
import hashlib, pickle, gzip

import numpy
import pandas

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec

from variation.variations import VariationsH5

import circadian
from introgressions_for_genes import calc_introgression_freqs_for_genes_from_variations
import introgressions
import passport
import pop_building
import genome
import rarefaction_vars
import snp_filtering


# upstream region 3kb annotated genes with most introgressions in Ec
# Solyc01g006680 AT3G20810 JMJD5/ JUMONJI DOMAIN CONTAINING 5 inferred from mutant phenotype: analysis of physiological response:Jones, et al. (2010)
# Solyc09g090100 AT1G04400 CRY2/ CRYPTOCHROME 2 inferred from expression pattern: Toth, et al. (2001)  
# Solyc10g086000 AT5G02810 PRR7/ PSEUDO-RESPONSE REGULATOR 7 traceable author statement: Matsushika, et al. (2000)  
# Solyc02g070530 AT4G18290 KAT2/ POTASSIUM CHANNEL IN ARABIDOPSIS THALIANA 2 inferred from mutant phenotype: analysis of physiological response: Lebaudy, et al. (2008)  
# Solyc06g051570 AT5G05660 NFXL2/ NFX1-LIKE 2 inferred from mutant phenotype: expression of another gene in a mutant background of this gene: Johansson, et al. (2011)  
# Solyc04g077970 AT1G27450 APT1/ ADENINE PHOSPHORIBOSYL TRANSFERASE 1 inferred from expression pattern:Transcript levels (e.g. RT-PCR):Li, et al. (2011)
# Solyc10g005080 AT2G46830 CCA1/ CIRCADIAN CLOCK ASSOCIATED 1 inferred from expression pattern: Transcript levels (e.g. RT-PCR): Li, et al. (2011)  

# upstream region 3kb mutant genes with most introgressions in Ec
# Solyc11g065120 AT5G63860 Ultraviolet-B receptor UVR8
# Solyc07g045480 AT5G35840 Phytochrome C 
# Solyc09g090100 AT1G04400 Cryptochrome-2


def get_introgression_freqs_for_genes(gene_ids, introgression_freqs_for_genes):
    gene_ids = list(gene_ids)
    return pandas.Series(numpy.array([introgression_freqs_for_genes.get(gene_id, 0) for gene_id in gene_ids]),
                         index=gene_ids)


def calc_diversities_around_genes(variations, genes, window_size=10000, allowed_missing_gts=0,
                                  ploidy=2, cache_dir=None):

    cache_path = None
    if cache_dir:
        key = ','.join(sorted(variations.samples))
        key += str(variations.num_variations)
        key += 'gene_ids' + str(list(genes.gene_ids))
        key += 'window_size' + str(window_size)
        key += 'allowed_missing_gts' + str(allowed_missing_gts)
        key += 'ploidy' + str(ploidy)
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('diversities_around genes.' + key + '.pickle')
        if cache_path.exists():
            return pickle.load(gzip.open(cache_path, 'rb'))

    res = rarefaction_vars.calc_pop_stats_per_var(variations, allowed_missing_gts=allowed_missing_gts,
                                                  ploidy=ploidy)
    exp_het_per_var = res['unbiased_exp_het']
    pos_index = variations.pos_index
    diversities = []
    gene_ids = []
    for gene_id in genes.gene_ids:
        gene = genes.get_gene(gene_id)
        mean_pos = int((gene['Start'] + gene['End']) / 2)
        start = int(mean_pos - (window_size / 2))
        if start < 0:
            start = 0
        end = int(mean_pos + (window_size / 2))
        chrom = gene['Chromosome'].encode()
        if chrom == b'SL4.0ch00':
            diversity = numpy.nan
        else:
            idx0 = pos_index.index_pos(chrom, start)
            idx1 = pos_index.index_pos(chrom, end) + 1
            diversity = numpy.nanmean(exp_het_per_var[idx0:idx1])
        diversities.append(diversity)
        gene_ids.append(gene_id)
    diversities = pandas.Series(diversities, index=gene_ids)

    if cache_dir:
        pickle.dump(diversities, gzip.open(cache_path, 'wb'))

    return diversities


def plot_introgressions_diversities(introgression_freqs_for_all_genes, introgression_freqs_for_circadian_genes, gene_diversities, fig,
                                    n_bins_introgressions=10,
                                    n_bins_diversities=10,
                                    introgression_range=(0, .9),
                                    diversity_range=(0, 0.4)):
    all_gene_diversities = gene_diversities.reindex(introgression_freqs_for_all_genes.index)

    spec = gridspec.GridSpec(ncols=3, nrows=3, figure=fig,
                             width_ratios=[0.8, 0.1, 0.1],
                             height_ratios=[0.8, 0.1, 0.1],
                             )
    axes = fig.add_subplot(spec[0, 0])
    axes.scatter(all_gene_diversities.values, introgression_freqs_for_all_genes, alpha=0.1)

    circadian_gene_diversities = gene_diversities.reindex(introgression_freqs_for_circadian_genes.index)
    axes.scatter(circadian_gene_diversities.values, introgression_freqs_for_circadian_genes, alpha=1)
    axes.set_xlim(*diversity_range)
    axes.set_ylim(*introgression_range)

    axes_introgression_all_hist = fig.add_subplot(spec[0, 1], sharey=axes)
    axes_introgression_all_hist.hist(introgression_freqs_for_all_genes.values, orientation="horizontal", bins=n_bins_introgressions, range=introgression_range)
    axes_introgression_all_hist.set_xlabel('all')

    axes_introgression_circadian_hist = fig.add_subplot(spec[0, 2], sharey=axes)
    axes_introgression_circadian_hist.hist(introgression_freqs_for_circadian_genes.values, orientation="horizontal", bins=n_bins_introgressions, range=introgression_range, color='blue')
    axes_introgression_circadian_hist.set_xlabel('circadian')

    axes_diversities_all_hist = fig.add_subplot(spec[1, 0], sharex=axes)
    axes_diversities_all_hist.hist(all_gene_diversities.values, bins=n_bins_diversities, range=diversity_range)
    axes_diversities_all_hist.set_ylabel('all')

    axes_diversities_circadian_hist = fig.add_subplot(spec[2, 0], sharex=axes)
    axes_diversities_circadian_hist.hist(circadian_gene_diversities.values, bins=n_bins_diversities, range=diversity_range, color='blue')
    axes_diversities_circadian_hist.set_ylabel('circadian')

    axes_diversities_circadian_hist.set_xlabel('diversities')
    axes.set_ylabel('introgression freq.')


def calc_distance_to_circadian_gene(genes, circadian_genes):
    all_genes_df = genes._genes_df
    circdian_genes_ids = list(circadian_genes.keys())
    circadian_genes_df = all_genes_df.reindex(circdian_genes_ids)

    all_genes_fhand = genes._write_bed(all_genes_df)
    circadian_genes_fhand = genes._write_bed(circadian_genes_df)

    cmd = ['bedtools', 'closest', '-b', circadian_genes_fhand.name, '-a', all_genes_fhand.name]
    process = subprocess.run(cmd, check=True, capture_output=True, text=True)

    dists = {}
    for line in process.stdout.splitlines():
        items = line.split()
        var_start = int(items[1])
        var_end = int(items[2])
        gene_start = int(items[5])
        gene_end = int(items[6])
        gene_id = items[3]

        dist_to_gene = min((abs(var_start - gene_start), abs(var_start - gene_end), abs(var_end - gene_start), abs(var_end - gene_end)))
        try:
            prev_dist = dists[gene_id]
            if prev_dist > dist_to_gene:
                dists[gene_id] = dist_to_gene
        except KeyError:
            dists[gene_id] = dist_to_gene
    dists = pandas.Series(dists)
    return dists


def plot_dists_to_circadian_genes(introgression_freqs_for_all_genes, gene_distances_to_circadian_gene, fig,
                                  n_bins=10, introgression_range=(0, .9)):
    dists_for_all_genes = gene_distances_to_circadian_gene.reindex(introgression_freqs_for_all_genes.index)

    spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig,
                             width_ratios=[0.85, 0.15],
                             height_ratios=[0.85, 0.15],
                             )
    axes = fig.add_subplot(spec[0, 0])
    axes.scatter(dists_for_all_genes.values, introgression_freqs_for_all_genes, alpha=0.1)

    axes_introgression_all_hist = fig.add_subplot(spec[0, 1], sharey=axes)
    axes_introgression_all_hist.hist(introgression_freqs_for_all_genes.values, orientation="horizontal", bins=n_bins, range=introgression_range)

    axes_dists_hist = fig.add_subplot(spec[1, 0], sharex=axes)
    axes_dists_hist.hist(dists_for_all_genes.values, bins=n_bins)


if __name__ == '__main__':

    debug = False
    cache_dir = config.CACHE_DIR

    diversity_range = 0, 0.1
    introgression_range = 0, 0.8
    n_bins_introgressions = 20
    n_bins_diversities = 10

    circadian_references = ['https://doi.org/10.1016/j.molp.2020.02.006',
                            'https://doi.org/10.1093/jxb/eru441',
                            'https://doi.org/10.3389/fpls.2020.00864',
                            'https://doi.org/10.1073/pnas.1801862115',
                            'https://doi.org/10.1073/pnas.1310631110',
                            'tair']
    reference_id = 'all_references'

    if True:
        circadian_references = ['https://doi.org/10.1016/j.molp.2020.02.006',
                                'https://doi.org/10.1093/jxb/eru441',
                                'https://doi.org/10.3389/fpls.2020.00864',
                                'https://doi.org/10.1073/pnas.1801862115',
                                'https://doi.org/10.1073/pnas.1310631110'
                                ]
        reference_id = 'papers'

    only_first = True

    if debug:
        vars_path = config.WORKING_H5
    else:
        vars_path = config.TIER1_H5_LOWQ_085

    variations = VariationsH5(str(vars_path), 'r')

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops_rank1 = pop_building.get_pops(pops_descriptions, passports)

    win_size = 10000
    method = 'mean_highest'

    genes = genome.Genes()

    founder_pop = 'slc_ma'

    target_and_introgression_source = [('slc_pe_n', config.RANK1, 'sp_pe'),
                                        ('slc_ec', config.RANK1, 'sp_ec'),
                                        ]

    only_first_text = 'only_first_hit' if only_first else 'all_blast_hists'

    out_dir = config.FIGURES_DIR
    out_dir.mkdir(exist_ok=True)

    pops_descriptions = {config.RANK1: ['slc_ma']}
    pops = pop_building.get_pops(pops_descriptions, passports)
    slc_ma_samples = pops['slc_ma']

    variations_for_slc_ma_diversities = snp_filtering.filter_variations(variations,
                                                                        samples_to_keep=slc_ma_samples,
                                                                        cache_dir=cache_dir)

    gene_diversities = calc_diversities_around_genes(variations_for_slc_ma_diversities,
                                                     window_size=500000,
                                                     allowed_missing_gts=20,
                                                     genes=genes,
                                                     cache_dir=cache_dir)

    tomato_genes = circadian.get_circadian_loci(only_first=only_first, only_references=circadian_references)

    if False:
        gene_distances_to_circadian_gene = calc_distance_to_circadian_gene(genes, tomato_genes)
        fig = Figure()
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)
        range_ = 0, 1
        axes.hist(gene_distances_to_circadian_gene, density=False)
        axes.legend()
        plot_path = config.TMP_DIR / f'dist_to_circadian_gene.{reference_id}.{only_first_text}.svg'
        fig.savefig(plot_path)

    for target_pop, target_pop_rank, introgression_source_pop in target_and_introgression_source:
        print(f'Introgressions in {target_pop}')

        samples_in_pop_with_possible_introgressions = pops_rank1[target_pop]

        samples_in_founder_pop = pops_rank1[founder_pop]
        samples_in_introgression_source_pop = pops_rank1[introgression_source_pop]

        introgession_config = {'samples_in_pop_with_introgressions': sorted(samples_in_pop_with_possible_introgressions),
                               'samples_in_founder_pop': sorted(samples_in_founder_pop),
                               'samples_in_introgression_source_pop': sorted(samples_in_introgression_source_pop),
                               'freq_threshold_to_consider_allele_present_in_founder_pop' : config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_PRESENT,
                               'freq_threshold_to_consider_allele_common_in_introgression_source_pop': config.FREQ_THRESHOLD_TO_CONSIDER_ALLELE_COMMON,
                              }
        introgression_freqs_for_genes = calc_introgression_freqs_for_genes_from_variations(variations,
                                                                                           introgession_config=introgession_config,
                                                                                           genes=genes,
                                                                                           upstream_region=config.UPSTREAM_REGION_IN_BP,
                                                                                           cache_dir=cache_dir
                                                                                           )

        introgression_freqs_for_circadian_genes = get_introgression_freqs_for_genes(tomato_genes.keys(),
                                                                                    introgression_freqs_for_genes)
        print(introgression_freqs_for_circadian_genes[introgression_freqs_for_circadian_genes > 0.4])

        introgression_freqs_for_all_genes = get_introgression_freqs_for_genes(genome.Genes().gene_ids,
                                                                              introgression_freqs_for_genes)

        fig = Figure((10, 10))
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        plot_introgressions_diversities(introgression_freqs_for_all_genes,
                                        introgression_freqs_for_circadian_genes,
                                        gene_diversities,
                                        fig,
                                        diversity_range=diversity_range,
                                        introgression_range=introgression_range,
                                        n_bins_diversities=n_bins_diversities,
                                        n_bins_introgressions=n_bins_introgressions,
                                        )
        plot_path = config.TMP_DIR / f'circadian_genes_introgressions_vs_diversities.{reference_id}.{only_first_text}.{target_pop}.svg'
        fig.savefig(plot_path)

        if False:
            fig = Figure((10, 10))
            FigureCanvas(fig) # Don't remove it or savefig will fail later
            plot_dists_to_circadian_genes(introgression_freqs_for_all_genes, gene_distances_to_circadian_gene, fig)
            plot_path = config.TMP_DIR / f'introgression_freqs_vs_dist_to_circadian_gene.{reference_id}.{only_first_text}.{target_pop}.svg'
            fig.savefig(plot_path)

        most_introgressed_genes = introgression_freqs_for_circadian_genes[introgression_freqs_for_circadian_genes > 0.1]
        for locus_id in most_introgressed_genes.index:
            print(locus_id, tomato_genes[locus_id].get('similar_to_arabidopsis', ''))

        fig = Figure()
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)
        range_ = 0, 1
        axes.hist(introgression_freqs_for_circadian_genes, label='circadian', density=True, range=range_, alpha=0.5)
        axes.hist(introgression_freqs_for_all_genes, label='all', density=True, range=range_, alpha=0.5)
        axes.legend()
        debug_text = '.debug' if debug else ''
        plot_path = config.TMP_DIR / f'circadian_genes_introgressions.{reference_id}.{only_first_text}{debug_text}.{target_pop}.svg'
        fig.savefig(plot_path)

        fig = Figure()
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)
        range_ = 0, 1
        axes.hist(introgression_freqs_for_circadian_genes, label='circadian', density=True, range=range_, alpha=0.5)
        axes.hist(introgression_freqs_for_all_genes, label='all', density=True, range=range_, alpha=0.5)
        axes.set_ylim((0, 0.5))
        axes.legend()
        plot_path = config.TMP_DIR / f'circadian_genes_introgressions.y_cut05.{reference_id}.{only_first_text}{debug_text}.{target_pop}.svg'
        fig.savefig(plot_path)
