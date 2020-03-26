import config

import re
from pandas import read_excel


if __name__ == '__main__':

    flanking_sequences_50bp_or_less_regex = '(.{,50})\[./.\](.{,50})'
    flanking_sequences_fasta = open(config.SOLCAP_FLANKING_SEQUENCES_FASTA, 'w')

    SolCAP_snps = read_excel(config.SOLCAP_ANOTATION_XLS)

    for id, flanking_sequence in zip(SolCAP_snps['SolCAP_SNP_ID'], SolCAP_snps['Flanking_Sequence']):
        flanking_sequences_50bp_or_less = re.search(flanking_sequences_50bp_or_less_regex, flanking_sequence)

        flanking_sequences_fasta.write('>' + id + '|left' + '\n')
        flanking_sequences_fasta.write(flanking_sequences_50bp_or_less.group(1) + '\n')
        flanking_sequences_fasta.write('>' + id + '|right' + '\n')
        flanking_sequences_fasta.write(flanking_sequences_50bp_or_less.group(2) + '\n')