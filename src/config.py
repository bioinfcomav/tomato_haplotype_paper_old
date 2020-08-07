
import sys
from pathlib import Path
import getpass

HOME_DIR = Path.home()
USER = getpass.getuser()

FINAL_VERSION = False

if FINAL_VERSION:
    SKIP_SAMPLES_WITH_NO_PASSPORT = False
    SKIP_SAMPLES_WITH_NO_GENOTYPE = False
else:
    SKIP_SAMPLES_WITH_NO_PASSPORT = True
    SKIP_SAMPLES_WITH_NO_GENOTYPE = True


if USER == 'jose':
    
    ANALYSES_DIR = HOME_DIR / 'magnet/analyses/'
    BASE_DIR = ANALYSES_DIR / 'ultimate_paper/'

    MODULE_DIR = BASE_DIR / 'src'
    sys.path.insert(0, str(MODULE_DIR))
    DEVEL_DIR = HOME_DIR / 'devel'
    VARIATION_DIR = DEVEL_DIR / 'variation5'
    sys.path.insert(0, str(VARIATION_DIR))

    TOMATO_GENOME_DIR = HOME_DIR / 'magnet/genomes/s_lycopersicum/sl40/'

    PYTHON2_BIN_FOR_FASTSTRUCTURE = ANALYSES_DIR / 'population_ultimate_env2/bin/python'
    FASTSTRUCTURE_BIN = ANALYSES_DIR / 'population_ultimate_env2/fastStructure-master/structure.py' 
    BEAGLE_JAR = HOME_DIR / 'soft/beagle/beagle.18May20.d20.jar'

elif USER == 'david':

    BASE_DIR = HOME_DIR / '/home/david/Disco/ultimate_paper/'
    MODULE_DIR = BASE_DIR / 'src'
    sys.path.insert(0, str(MODULE_DIR))

    DEVEL_DIR = HOME_DIR / '/home/david/Development'
    VARIATION_DIR = DEVEL_DIR / 'variation5'
    sys.path.insert(0, str(VARIATION_DIR))

    TOMATO_GENOME_DIR = HOME_DIR / 'Disco/genomes/s_lycopersicum/sl40/'

elif USER == 'jope':

    BASE_DIR = HOME_DIR / '/home/jope/tomato/ultimate_paper'
    MODULE_DIR = BASE_DIR / 'src'
    sys.path.insert(0, str(MODULE_DIR))

    DEVEL_DIR = HOME_DIR / 'devel3'
    VARIATION_DIR = DEVEL_DIR / 'variation5'
    sys.path.insert(0, str(VARIATION_DIR))

    TOMATO_GENOME_DIR = HOME_DIR / 'genomes/tomato/SL4.00/'
    TOMATO_GENOME_FASTA = HOME_DIR / 'S_lycopersicum_chromosomes.4.00.fa'

TOMATO_GENOME_FASTA = HOME_DIR / 'S_lycopersicum_chromosomes.4.00.fa'

SOURCE_CODE_DIR = BASE_DIR / 'src'
TREE_MIX_PLOTTING_R_SOURCE = SOURCE_CODE_DIR / 'plotting_funcs.R'

SOURCE_DATA_DIR = BASE_DIR / 'source_data'

PASSPORTS_CSV = SOURCE_DATA_DIR / 'passports.csv'
WORKING_PASSPORTS_CSV = SOURCE_DATA_DIR / 'working_passports.csv'

MAP_BACKGROUNDS_DIR = SOURCE_DATA_DIR / 'map_backgrounds'

NE_BACKGROUND_TIFF = MAP_BACKGROUNDS_DIR / 'NE1_50M_SR_W' / 'NE1_50M_SR_W.tif'
NE_BACKGROUND_CUT_PNG = MAP_BACKGROUNDS_DIR / 'ne1_50m_sr.cut.png'

FIGURES_DIR = BASE_DIR / 'figures'

GEOGRAPHIC_FIGURE_DIR = FIGURES_DIR / 'geographic_map'

HYPOTHESIS_PNG = GEOGRAPHIC_FIGURE_DIR / 'hypothesis.png'

CACHE_DIR = BASE_DIR / 'caches'
CACHE_DIR.mkdir(exist_ok=True)

SNPS_DIR = BASE_DIR / 'snps'
FILTERING_STATS_DIR = SNPS_DIR / 'filtering_stats'
IMPUTATION_DIR = BASE_DIR / 'imputation'

ORIG_H5 = SNPS_DIR / 'tomato_genomic-20200204.h5'
TIER1_H5_LOWQ_085 = SNPS_DIR / 'tomato_genomic-20200204.tier1.low_qual_085.h5'
#TIER1_H5_LOWQ_090 = SNPS_DIR / 'tomato_genomic-20200204.tier1.low_qual_090.h5'
WORKING_H5 = SNPS_DIR / 'tomato_genomic-20200204.working_set.h5'
WORKING_VCF = SNPS_DIR / 'tomato_genomic-20200204.working_set.vcf'
WORKING_PHASED_H5 = SNPS_DIR / 'tomato_genomic-20200204.working_set.phased.h5'
WORKING_PHASED_AND_IMPUTED_H5 = SNPS_DIR / 'tomato_genomic-20200204.working_set.phased_and_imputed.h5'

#LOW_QUAL_SAMPLES_090 = SOURCE_DATA_DIR / 'low_quality_samples_090.txt'
LOW_QUAL_SAMPLES_085 = SOURCE_DATA_DIR / 'low_quality_samples_085.txt'

LOW_QUAL_MIN_SAMPLE_CALLED_RATE = 0.85
LOW_QUAL_N_BINS = 60

TIER1_MIN_GT_DP = 5
TIER1_MIN_SNP_QUAL = 50
TIER1_MAX_MAC = 2
TIER1_MIN_CALLED = 0.60
TIER1_MAX_HET = 0.10
TIER1_MAX_HET_MIN_CALL_DP = 20
TIER1_MAX_MAF = None

TIER2_MAX_MAC = 3
TIER2_MIN_CALLED = 0.90
TIER2_MAX_HET = 0.10
TIER2_MAX_HET_MIN_CALL_DP = 20
TIER2_MAX_MAF = None
TIER2_PCA_MAX_MAF = 0.99

SOLCAP_SOURCE_DIR = SOURCE_DATA_DIR / 'solcap'
SOLCAP_ANOTATION_XLS = SOLCAP_SOURCE_DIR / 'tomato_ng_snp_10k_infinium_annotation_v2.10_FINAL_w_solcap_ids.xls'
SOLCAP_GENETIC_MAP = SOLCAP_SOURCE_DIR / 'SolCap_genetic_map.csv'
BEAGLE_MAP = SOLCAP_SOURCE_DIR / 'snps_interpolated.beagle.map'

RELEVANT_FIELDS = ['/calls/DP', '/calls/GT', '/variations/chrom',
                   '/variations/pos', '/variations/ref',
                   '/variations/alt']

CLASSIFICATIONS = BASE_DIR / 'classifications.csv'
CLASSIFICATIONS_MORPHOLOGICAL = BASE_DIR / 'classifications_morphological.csv'
WORKING_CLASSIFICATION_CSV = BASE_DIR / 'working_classifications.csv'
CLASSIFICATION_HISTORY_DIR = BASE_DIR / 'old_classifications'

CLASSIFICATION_RANKS = ['rank1', 'rank2', 'morpho_type']

RANK1 = 'classification', 'rank1'
RANK2 = 'classification', 'rank2'
RAZIFARD = 'razifard_classification',
KEEP = 'keep'
REMOVE = 'remove'
ALL_POPS = 'all_pops'

SNPS_DIR = BASE_DIR / 'variations'

TIER1_PHASED_LOW_QUAL_09_MISSING_085 = SNPS_DIR / \
    '20190705.tier1.euchromatic.missing085_low_qual_09.phased.h5'
TIER1_PHASED_AND_IMPUTED_LOW_QUAL_09_MISSING_085 = SNPS_DIR / \
    '20190705.tier1.euchromatic.missing085_low_qual_09.phased_and_imputed.h5'

MULTIVAR_DIR = BASE_DIR / 'multivar'

HAPLO_PCOA_DIR = BASE_DIR / 'haplo_pcoa'

MIN_NUM_SNPS_FOR_HAPLO_IN_PCA = 20
HAPLO_WIN_SIZE = 5e5

SOLCAP_DIR = BASE_DIR / 'solcap'

RAREFACTION_DIR = BASE_DIR / 'rarefaction'
RAREFACTION_VARS_DIR = RAREFACTION_DIR / 'vars'
RAREFACTION_HAPLO_DIR = RAREFACTION_DIR / 'haplo'

CLASSIFICATION_CONFIG = {'thinning_dist_threshold': 0.00030,
                         'method': 'agglomerative',
                         'n_clusters': 3}
CLASSIFICATION_OUTLIER_CONFIG = {'method': 'elliptic_envelope',
                                 'contamination': 0.2}
SUPERVISED_CLASSIFICATION_CONFIG = {'prob_threshold': 0.99,
                                    'classifier': 'kneighbors',
                                    'n_neighbors': 30}
CLASSIFICATION_REFERENCES = {'SL4.0ch01%610462%ts-554%1': 'sl',
                             'SL4.0ch01%610462%ts-450%1': 'sp_peru',
                             'SL4.0ch01%610462%bgv007339%1': 'sp_ecu'}
OUTLIER_CONFIGS = [{'method': 'isolation_forest', 'contamination': 0.070,
                    'thinning_dist_threshold': 0.0015}]
N_DIMS_TO_KEEP = 3

DIVERSITIES_DIR = BASE_DIR / 'diversities'
DIVERSITIES_VAR_DIR = DIVERSITIES_DIR / 'vars'
DIVERSITIES_HAPLO_DIR = DIVERSITIES_DIR / 'haplo'

TREE_MIX_DIR = BASE_DIR / 'tree_mix'

FASTSTRUCTURE_DIR = BASE_DIR / 'faststructure'
FASTSTRUCTURE_RUN_DIR = FASTSTRUCTURE_DIR / 'runs'
FASTSTRUCTURE_RESULT_BASE_FNAME = 'structure_out'
FASTSTRUCTURE_PLINK_DIR = FASTSTRUCTURE_DIR / 'plink_files'
FASTSTRUCTURE_PLOT_DIR = FASTSTRUCTURE_DIR / 'plots'
FASTSTRUCTURE_VS_HAPLO_PLOT_DIR = FASTSTRUCTURE_PLOT_DIR / 'faststructure_vs_haplo_classes'

PLINK_BIN = 'plink1.9'

LD_DIR = BASE_DIR / 'ld'

TMP_DIR = BASE_DIR / 'tmp'

FREQ_THRESHOLD_TO_CONSIDER_ALLELE_PRESENT = 0.01
FREQ_THRESHOLD_TO_CONSIDER_ALLELE_COMMON = 0.1

THERE_AND_BACK_DIR = BASE_DIR / 'there_and_back'

CDNA_FASTA = TOMATO_GENOME_DIR / 'ITAG4.1_cDNA.fasta '
GENE_MODELS_GFF = TOMATO_GENOME_DIR / 'ITAG4.1_gene_models.gff.gz'

MORPHOLOGICAL_SOURCE_DATA = SOURCE_DATA_DIR / 'morphological_characterization.xlsx'

MORPHOLOGICAL_DIR = BASE_DIR / 'morphology'

SELECTIVE_SWEEP_DIR = BASE_DIR / 'selective_sweep'

RELEVANT_GENES_DIR = BASE_DIR / 'relevant_genes'

BAM_LIST = SOURCE_DATA_DIR / 'bam.list'

BAM_STATS_DIR = SOURCE_DATA_DIR / 'bamstats'

RAZIFARD_CLASSIFICATION = SOURCE_DATA_DIR / 'razifard_classification.csv'

ECOSYSTEM_DATA = SOURCE_DATA_DIR / 'ecosystem_data.csv'

TABLE_MAPPING_STATS = FIGURES_DIR / 'suppl_table_1_mapping_stats.csv'
TABLE_PASSPORT_AND_MORPHOLOGICAL_DATA = FIGURES_DIR / 'suppl_table_2_passport_and_morphological_data.csv'
TABLE_GENES_WITH_MANY_INTROGRESSIONS = FIGURES_DIR / 'suppl_table_3_genes_with_many_introgressions.csv'

FIG_HAPLO_PCA = FIGURES_DIR / 'fig_1_haplo_pca.svg'
FIG_ACC_PCA = FIGURES_DIR / 'fig_2_acc_pca.svg'
FIG_HAPLO_POP_COMPOSITION = FIGURES_DIR / 'fig_3_haplo_pop_composition.png'
FIG_DIVERSITIES_BY_HAPLO_TYPE = FIGURES_DIR / 'fig_4_diversities_by_haplo_kind.svg'
FIG_DIVERSITIES_AND_LD = FIGURES_DIR / 'fig_5_diversities_and_ld.svg'
FIG_THERE_AND_BACK = FIGURES_DIR / 'fig_6_slc_ma_genomic_diversity_vs_andean_slc_introgressions.svg'
FIG_HSL_SHARED_HAPLOS = FIGURES_DIR / 'fig_7_shared_hsl_haplos.svg'
FIG_MORPHOLOGICAL_ANALYSIS = FIGURES_DIR / 'fig_8_morphological_analaysis.svg'

FIG_NUM_HAPLO_TYPES = FIGURES_DIR / 'suppl_fig_1_num_haplo_types.svg'
FIG_HAPLO_PCA_PER_SAMPLE = FIGURES_DIR / 'suppl_fig_2_haplo_pca_per_sample.png'
FIG_ACC_HAPLO_FREQS_AND_STRUCTURE = FIGURES_DIR / 'suppl_fig_3_acc_haplo_freqs_and_structure.svg'
FIG_ACC_PCA_HIERALCHICAL = FIGURES_DIR / 'suppl_fig_4_hieralchical_acc_pcas.svg'
FIG_RAZIFARD_VS_CURRENT_CLASSIFICATION = FIGURES_DIR / 'suppl_fig_5_suppl_razifard_classification_comparison.svg'
FIG_RAREFACTION = FIGURES_DIR / 'suppl_fig_7_rarefaction.svg'
FIG_PI_DISTRIBUTIONS = FIGURES_DIR / 'suppl_fig_8_pi_distributions.svg'
FIG_HSL_SHARED_HAPLOS_2 = FIGURES_DIR / 'suppl_fig_9_shared_hsl_haplos.svg'
FIG_HAPLO_FREQS = FIGURES_DIR / 'suppl_fig_10_haplo_freqs.svg'
FIG_MORPHOLOGICAL_TYPES_CHARACTERIZATION = FIGURES_DIR / 'suppl_fig_11_morphological_types_characterization.svg'
FIG_COLLECTING_SOURCES = FIGURES_DIR / 'suppl_fig_12_collecting_sources.svg'
