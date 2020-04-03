
import config

import re
import itertools
import csv
import hashlib, pickle
from collections import defaultdict
from pprint import pprint

import numpy
import pandas
from scipy import interpolate

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.patches as patches

import locate_seq_using_blast

IUPAC_SNP_CODES = {'AG': 'R',
                   'GA': 'R',
                   'CT': 'Y',
                   'TC': 'Y',
                   'CG': 'S',
                   'GC': 'S',
                   'AT': 'W',
                   'TA': 'W',
                   'GT': 'K',
                   'TG': 'K',
                   'AC': 'M',
                   'CA': 'M',
                   'CGT': 'B',
                   'CTG': 'B',
                   'GCT': 'B',
                   'GTC': 'B',
                   'TGC': 'B',
                   'TCG': 'B',
                   'AGT': 'D',
                   'ATG': 'D',
                   'TAG': 'D',
                   'TGA': 'D',
                   'GTA': 'D',
                   'GAT': 'D',
                   'ACT': 'H',
                   'ATC': 'H',
                   'TAC': 'H',
                   'TCA': 'H',
                   'CTA': 'H',
                   'CAT': 'H',
                   'ACG': 'V',
                   'AGC': 'V',
                   'CAG': 'V',
                   'CGA': 'V',
                   'GCA': 'V',
                   'GAC': 'V',
                   }


NUM_LINEAL_SEGMENTS = {'SL4.0ch01': 3,
                       'SL4.0ch02': 3,
                       'SL4.0ch03': 5,
                       'SL4.0ch04': 5,
                       'SL4.0ch05': 5,
                       'SL4.0ch06': 4,
                       'SL4.0ch07': 4,
                       'SL4.0ch08': 4,
                       'SL4.0ch09': 5,
                       'SL4.0ch10': 3,
                       'SL4.0ch11': 6,
                       'SL4.0ch12': 5}


def get_IUPAC_snp_code(*args):
    if not 2 <= len(args) <= 3:
        raise ValueError('only acepts 2 or 3 snp variants')

    snp_variants = list(args)
    #snp_variants.sort()
    snp_variants = ''.join(snp_variants)

    return IUPAC_SNP_CODES[snp_variants]


def _replace_snp_nucleotides(seq_with_snp):
    flanking_sequence_with_variants_regex = '(.+)\[(.)/(.)\](.+)'
    flanking_sequences_and_snp_variants = re.search(flanking_sequence_with_variants_regex, seq_with_snp)
    snp = get_IUPAC_snp_code(flanking_sequences_and_snp_variants.group(2), flanking_sequences_and_snp_variants.group(3))
    flanking_sequence = flanking_sequences_and_snp_variants.group(1) + snp + flanking_sequences_and_snp_variants.group(4)

    return flanking_sequence


def read_solcap_genetic_map():
    markers = {}
    for row in csv.DictReader(config.SOLCAP_GENETIC_MAP.open('rt')):
        markers[row['snp']] = {'chrom': row['gen_chr'], 'genet_loc': row['gen_pos']}
    return markers


def get_solcap_markers(approx_phys_loc=False, cache_dir=None):

    if cache_dir:
        key = str(approx_phys_loc)
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('solcap_markers' + key + '.pickle')
        if cache_path.exists():
            return pickle.load(cache_path.open('rb'))

    if not approx_phys_loc:
        raise NotImplementedError()

    genetic_map = read_solcap_genetic_map()

    solcap_df = pandas.read_excel(config.SOLCAP_ANOTATION_XLS)
    markers = {}
    for _, row in solcap_df.iterrows():
        marker_id = row['SolCAP_SNP_ID']

        chrom25 = row['Chromosome3']

        try:
            genet_loc = genetic_map[marker_id]
        except KeyError:
            continue

        seq = _replace_snp_nucleotides(row['Flanking_Sequence'])
        try:
            locations = locate_seq_using_blast.locate_sequence(seq, config.TOMATO_GENOME_FASTA, check_uniq=True)
        except RuntimeError:
            continue
        if not locations:
            continue

        location = locations[0]
        chrom = location['subject']

        if chrom25 != 'unknown' and int(chrom25[-2:]) != int(chrom[-2:]):
            continue
        
        location = int((location['end'] + location['start']) / 2)

        if int(genet_loc['chrom']) != int(chrom[-2:]):
            continue
        genet_loc = float(genet_loc['genet_loc'])

        markers[marker_id] = {'chrom': chrom,
                              'phys_loc': location,
                              'genet_loc': genet_loc
                             }
    if cache_dir:
        pickle.dump(markers, cache_path.open('wb'))

    return markers


def collect_locs_per_chrom(markers):

    locs = defaultdict(list)
    for marker in markers.values():
        locs[marker['chrom']].append((marker['phys_loc'], marker['genet_loc']))
    
    for chrom in locs.keys():
        locs[chrom].sort(key=lambda x: x[0])

    genet_locs = {}
    phys_locs = {}
    for chrom in locs.keys():
        genet_locs[chrom] = [loc[1] for loc in locs[chrom]]
        phys_locs[chrom] = [loc[0] for loc in locs[chrom]]

    return genet_locs, phys_locs


def calc_recomb_rate(phys_locs, model, phys_win_size):
    phys_locs = numpy.array(phys_locs)
    phys_win_size = phys_win_size / 2

    phys_locs_minus_h = phys_locs - phys_win_size
    phys_locs_plus_h = phys_locs + phys_win_size

    # we don't want to extrapolate
    interpolable_phys_locs = numpy.logical_and(phys_locs_minus_h >= numpy.min(phys_locs),
                                               phys_locs_plus_h <= numpy.max(phys_locs))
    phys_locs_minus_h = phys_locs_minus_h[interpolable_phys_locs]
    phys_locs_plus_h = phys_locs_plus_h[interpolable_phys_locs]
    phys_locs = phys_locs[interpolable_phys_locs]

    recomb_rate = (model(phys_locs_plus_h) - model(phys_locs_minus_h)) / (2 * phys_win_size)

    return phys_locs, recomb_rate


def plot_genet_vs_phys_loc(markers, out_dir, models=None, euchromatic_regions=None):

    if models is None:
        models = {}

    if euchromatic_regions is None:
        euchromatic_regions = {}

    genet_dists, phys_locs = collect_locs_per_chrom(markers)

    chroms = sorted(genet_dists.keys())
    for chrom in chroms:
        fig = Figure()
        FigureCanvas(fig) # Don't remove it or savefig will fail later

        axes = fig.add_subplot(311)
        recomb_rate_axes = fig.add_subplot(312, sharex=axes)
        modeled_axes = fig.add_subplot(313, sharex=axes)
        x_ticks_kargs = {'axis': 'x', 'labelbottom': False, 'labelright': False,
                         'bottom': False}
        axes.tick_params(**x_ticks_kargs)
        recomb_rate_axes.tick_params(**x_ticks_kargs)
        modeled_axes.set_xlabel('Chrom. location (bp)')
        modeled_axes.set_ylabel('Genet. dist.')
        axes.set_ylabel('Genet. dist.')
        recomb_rate_axes.set_ylabel('Recomb. rate (per bp)')

        model = models.get(chrom)

        y_values = genet_dists[chrom]

        axes.scatter(phys_locs[chrom], y_values, zorder=10)
        
        x_values = numpy.linspace(numpy.min(phys_locs[chrom]),
                                    numpy.max(phys_locs[chrom]),
                                    1000)
        y_values = model(x_values)

        modeled_axes.plot(x_values, y_values, zorder=20)
    
        x_values = numpy.linspace(numpy.min(phys_locs[chrom]),
                                        numpy.max(phys_locs[chrom]),
                                        1000)
        x_values, y_values = calc_recomb_rate(x_values, model, 10000)
        recomb_rate_axes.plot(x_values, y_values, zorder=20, label='recomb_rate')

        chrom_regions = euchromatic_regions.get(chrom)

        if chrom_regions:
            for region in chrom_regions:
                width = region[1] - region[0]
                color = '#aaaaaa' if region[2] else '#cccccc'
                for axes_ in [axes, recomb_rate_axes, modeled_axes]:
                    height = axes_.get_ylim()[1] - axes_.get_ylim()[0]
                    rect = patches.Rectangle((region[0], 0), width, height, edgecolor=None, facecolor=color)
                    axes_.add_patch(rect)

        plot_path = out_dir / f'genet_vs_phys_{chrom}.svg'
        fig.tight_layout()
        fig.savefig(str(plot_path))


def demand_increase_is_monotonic(vector):
    prev_values = numpy.array([0] + list(vector[:-1]))
    increment = vector - prev_values
    increment[increment < 0] = 0
    monotonic_increased_vector = numpy.cumsum(increment)
    return monotonic_increased_vector


class MonotonicIncrease():
    def __init__(self, funct):
        self.funct = funct

    def __call__(self, *args, **kwargs):
        return demand_increase_is_monotonic(self.funct(*args, **kwargs))


def fit_markers(markers, k=3, s=100, der=0):
    genet_locs, phys_locs = collect_locs_per_chrom(markers)
    chroms = sorted(genet_locs.keys())

    models = {}
    for chrom in chroms:
        chrom_genet_dists = genet_locs[chrom]
        chrom_phys_locs = phys_locs[chrom]

        tck = interpolate.splrep(chrom_phys_locs, chrom_genet_dists, k=k, s=s)
        chrom_genet_dists = interpolate.splev(chrom_phys_locs, tck, der=der)

        chrom_genet_dists = demand_increase_is_monotonic(chrom_genet_dists)

        fitted_funct = MonotonicIncrease(interpolate.interp1d(chrom_phys_locs, chrom_genet_dists, kind='linear'))
        models[chrom] = fitted_funct
    return models


def determine_eucrohomatic_regions(markers, models, win_size, recomb_rate_threshold):

    genet_locs, phys_locs = collect_locs_per_chrom(markers)
    chroms = sorted(genet_locs.keys())
    euchromatic_regions = {}
    for chrom in chroms:

        x_values, recomb_rate = calc_recomb_rate(phys_locs[chrom],
                                                 models[chrom],
                                                 win_size)

        are_euchromatic = recomb_rate >= recomb_rate_threshold
        current_region_start = None
        current_state = None
        last_x_value = None
        regions = []
        for x_value, is_euchromatic in zip(x_values, are_euchromatic):
            if current_region_start is None:
                current_region_start = x_value
                current_state = is_euchromatic
                last_x_value = None
                continue
            if last_x_value and is_euchromatic != current_state:
                regions.append((current_region_start, last_x_value, current_state))
                current_region_start = x_value
                current_state = is_euchromatic

            last_x_value = x_value
        regions.append((current_region_start, last_x_value, current_state))
        euchromatic_regions[chrom] = regions
    return euchromatic_regions


def generate_interpolated_map():
    markers = get_solcap_markers(approx_phys_loc=True,
                                 cache_dir=config.CACHE_DIR)
    models = fit_markers(markers)

    genet_locs, phys_locs = collect_locs_per_chrom(markers)

    locs = defaultdict(list)
    for chrom, model in models.items():
        phys_locs_to_interpolate = numpy.linspace(numpy.min(phys_locs[chrom]),
                                                  numpy.max(phys_locs[chrom]),
                                                  1000)

        interpolated_genet_dists = model(phys_locs_to_interpolate)
    
        locs[chrom] = list(zip(phys_locs_to_interpolate, interpolated_genet_dists))

    return locs


if __name__ == '__main__':

    generate_interpolated_map()

    markers = get_solcap_markers(approx_phys_loc=True,
                                 cache_dir=config.CACHE_DIR)

    models = fit_markers(markers)

    euchromatic_regions = determine_eucrohomatic_regions(markers, models,
                                                         win_size=1000,
                                                         recomb_rate_threshold=1e-6)
    print('chromatic regions')
    print([(chrom, start, end) for chrom, regions in euchromatic_regions.items() for start, end, is_euchromatic in regions if is_euchromatic])

    out_dir = config.SOLCAP_DIR
    out_dir.mkdir(exist_ok=True)
    out_dir /= 'chroms'
    out_dir.mkdir(exist_ok=True)

    plot_genet_vs_phys_loc(markers, out_dir, models, euchromatic_regions)