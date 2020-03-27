import config

import re
from pandas import read_excel


IUPAC_SNP_CODES = {'AG': 'R',
                   'CT': 'Y',
                   'CG': 'S',
                   'AT': 'W',
                   'GT': 'K',
                   'AC': 'M',
                   'CGT': 'B',
                   'AGT': 'D',
                   'ACT': 'H',
                   'ACG': 'V'}


def get_IUPAC_snp_code(*args):
    if not 2 <= len(args) <= 3:
        raise ValueError('only acepts 2 or 3 snp variants')

    snp_variants = list(args)
    snp_variants.sort()
    snp_variants = ''.join(snp_variants)

    return IUPAC_SNP_CODES[snp_variants]


def export_flanking_sequenes_with_random_var_to_fasta(ids, sequences, out_fhand):
    flanking_sequence_with_random_var_regex = '(.+)\[(.)/(.)\](.+)'

    for id_, flanking_sequence in zip(ids, sequences):
        flanking_sequences_and_snp_variants = re.search(flanking_sequence_with_random_var_regex, flanking_sequence)
        snp = get_IUPAC_snp_code(flanking_sequences_and_snp_variants.group(2), flanking_sequences_and_snp_variants.group(3))
        flanking_sequence = flanking_sequences_and_snp_variants.group(1) + snp + flanking_sequences_and_snp_variants.group(4)

        out_fhand.write('>' + id_ + '\n')
        out_fhand.write(flanking_sequence + '\n')


def export_left_and_right_flanking_sequences_to_fasta(ids, sequences, out_fhand):
    flanking_sequences_50bp_or_less_regex = '(.{,50})\[./.\](.{,50})'

    for id_, flanking_sequence in zip(ids, sequences):
        flanking_sequences_50bp_or_less = re.search(flanking_sequences_50bp_or_less_regex, flanking_sequence)

        out_fhand.write('>' + id_ + '|left' + '\n')
        out_fhand.write(flanking_sequences_50bp_or_less.group(1) + '\n')
        out_fhand.write('>' + id_ + '|right' + '\n')
        out_fhand.write(flanking_sequences_50bp_or_less.group(2) + '\n')


if __name__ == '__main__':

    SolCAP_snps = read_excel(config.SOLCAP_ANOTATION_XLS)
    snp_ids = SolCAP_snps['SolCAP_SNP_ID']
    flanking_secuences = SolCAP_snps['Flanking_Sequence']
    flanking_sequences_fasta = open(config.SOLCAP_FLANKING_SEQUENCES_FASTA, 'w')

    export_flanking_sequenes_with_random_var_to_fasta(snp_ids, flanking_secuences,
                                                       flanking_sequences_fasta)