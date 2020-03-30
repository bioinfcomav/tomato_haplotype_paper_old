
import config

import re
import itertools
import csv

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


def get_solcap_markers():

    genetic_map = read_solcap_genetic_map()

    solcap_df = pandas.read_excel(config.SOLCAP_ANOTATION_XLS)
    markers = {}
    for _, row in solcap_df.iterrows():
        marker_id = row['SolCAP_SNP_ID']

        chrom25 = row['Chromosome3']

        seq = _replace_snp_nucleotides(row['Flanking_Sequence'])
        try:
            locations = locate_seq_using_blast.locate_sequence(seq, config.TOMATO_GENOME_FASTA, check_uniq=True)
        except RuntimeError:
            continue
        if not locations:
            continue

        location = locations[0]
        chrom = location['subject']

        if chrom25[-2:] != chrom[-2:]:
            continue
        
        location = int((location['end'] + location['start']) / 2)

        try:
            genet_loc = genetic_map[marker_id]
        except KeyError:
            continue

        if genet_loc['chrom'] != chrom[-2:]:
            continue
        genet_loc = genet_loc['genet_loc']
        print(chrom, location, genet_loc)

        markers[marker_id] = {'chrom': chrom,
                              'phys_loc': location,
                              'genet_loc': genet_loc
                             }
    return markers


if __name__ == '__main__':
    get_solcap_markers()