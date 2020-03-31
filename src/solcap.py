
import config

import re
import itertools
import csv
import hashlib, pickle
from collections import defaultdict

import numpy
import pandas

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

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


def _collect_locs_per_chrom(markers):
    genet_locs = defaultdict(list)
    phys_locs = defaultdict(list)
    for marker in markers.values():
        genet_locs[marker['chrom']].append(marker['genet_loc'])
        phys_locs[marker['chrom']].append(marker['phys_loc'])
    return genet_locs, phys_locs


def plot_genet_vs_phys_loc(markers, out_dir, fitted_phys_dists=None, fitted_genet_dists=None):

    if fitted_genet_dists is None:
        fitted_genet_dists = {}
    if fitted_phys_dists is None:
        fitted_phys_dists = {}

    genet_locs, phys_locs = _collect_locs_per_chrom(markers)

    chroms = sorted(genet_locs.keys())
    for chrom in chroms:
        fig = Figure()
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111)

        axes.scatter(phys_locs[chrom], genet_locs[chrom])

        fitted_x_vals = fitted_phys_dists.get(chrom)
        fitted_y_vals = fitted_genet_dists.get(chrom)
        if fitted_x_vals is not None:
            axes.plot(fitted_x_vals, fitted_y_vals)

        plot_path = out_dir / f'genet_vs_phys_{chrom}.svg'
        fig.tight_layout()
        fig.savefig(str(plot_path))


def segmented_lineal_interpolation3(x_values,
                                    xe0, xe1, xe2, xe3,
                                    ye0, ye1, ye2, ye3):
    x_model_edges = [xe0, xe1, xe2, xe3]
    y_model_edges = [ye0, ye1, ye2, ye3]
    return segmented_lineal_interpolation(x_values, x_model_edges, y_model_edges)


def segmented_lineal_interpolation(x_values, x_model_edges, y_model_edges):
    y_values = []
    
    mask = x_values < x_model_edges[0]
    if numpy.any(mask):
        y_values = [y_model_edges[0]] * numpy.sum(mask)

    first = True
    for (x0, x1), (y0, y1) in zip(zip(x_model_edges[:-1], x_model_edges[1:]), zip(y_model_edges[:-1], y_model_edges[1:])):

        if first:
            mask = numpy.logical_and(x_values >= x0, x_values <= x1)
            first = False
        else:
            mask = numpy.logical_and(x_values > x0, x_values <= x1)
        segment_x_values = x_values[mask]
        slope = (y1 - y0) / (x1 - x0)
        orig = y0 - x0 * slope
        segment_y_values = (y0 - slope * x0) + (slope * segment_x_values)
        y_values.extend(segment_y_values)

    mask = x_values > x_model_edges[-1]
    if numpy.any(mask):
        y_values.extend([y_model_edges[-1]] * numpy.sum(mask))

    y_values = numpy.array(y_values)
    return y_values


def calc_genet_dists_from_model(model, phys_dists):
    chroms = sorted(model.keys())

    genet_dists = {}
    for chrom in chroms:
        x_model_edges = model[chrom][0]
        y_model_edges = model[chrom][1]
        x_values = phys_dists[chrom]
        y_values = segmented_lineal_interpolation(x_values, x_model_edges, y_model_edges)
        genet_dists[chrom] = y_values
    return genet_dists


if __name__ == '__main__':
    markers = get_solcap_markers(approx_phys_loc=True,
                                 cache_dir=config.CACHE_DIR)
    print(len(markers))
    print(sorted({marker['chrom'] for marker in markers.values()}))

    model = {'SL4.0ch01': ([0, 0.8e7, 7e7, 8.8e7], [0, 42, 42, 130])}
    phys_dists = {'SL4.0ch01': numpy.linspace(0, 8.8e7, 100)}
    genet_dists = calc_genet_dists_from_model(model, phys_dists)

    out_dir = config.SOLCAP_DIR
    out_dir.mkdir(exist_ok=True)
    out_dir /= 'chroms'
    out_dir.mkdir(exist_ok=True)

    plot_genet_vs_phys_loc(markers, out_dir, phys_dists, genet_dists)
