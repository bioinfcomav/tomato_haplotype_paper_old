
import config

import re
import itertools
import csv
import hashlib, pickle

import pandas

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
        genet_loc = genet_loc['genet_loc']

        markers[marker_id] = {'chrom': chrom,
                              'phys_loc': location,
                              'genet_loc': genet_loc
                             }
    if cache_dir:
        pickle.dump(markers, cache_path.open('wb'))

    return markers


if __name__ == '__main__':
    markers = get_solcap_markers(approx_phys_loc=True,
                                 cache_dir=config.CACHE_DIR)
    print(len(markers))
    print({marker['chrom'] for marker in markers.values()})

    #model = calculate_genet_loc_model(markers)
    #calc_genet_dists_from_model(model, phys_dists{'chrom01': phys_dists})

    #plot_genet_vs_phys_loc(markers, plot_path, model=None)

    