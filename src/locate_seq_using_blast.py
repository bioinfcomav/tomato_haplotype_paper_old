
import config

from pprint import pprint

import blast
import ranges


def locate_sequence(seq, genome, check_uniq=True, evalue_threshold=1e-40,
                   hsp_align_len_threshold=None, hsp_identity_threshold=None):

    blast_dir = config.CACHE_DIR / 'blast'
    blast_dir.mkdir(exist_ok=True)

    res = blast.prepare_blast_db(genome, blast_dir, db_type='nucl',
                                 skip_if_exists=True)
    db_path = res['db_path']

    res = blast.blast_seqs([{'seq': seq, 'name': 'seq'}],
                            db_path, blast_program='blastn',
                            evalue_threshold=evalue_threshold)

    subjects = list(res['seq'].keys())
    if check_uniq:
        if len(subjects) > 1:
            raise RuntimeError('Sequence found in several subjects')

    locations = []
    for subject, hsps in res['seq'].items():
        hsps = blast.filter_hsps_by_align_len(hsps, hsp_align_len_threshold)
        hsps = blast.filter_hsps_by_identity(hsps, hsp_identity_threshold)

        subject_ranges = []
        for hsp in hsps:
            subject_ranges.append({'chrom': subject,
                                   'start': hsp['subject_start'],
                                   'end': hsp['subject_end'],
                                   'strand': '+' if hsp['subject_strand'] == 1 else '-'})
        subject_ranges = ranges.merge_ranges(subject_ranges)
        for range_ in subject_ranges:
            locations.append({'subject': subject,
                              'start': range_['start'],
                              'end': range_['end'],
                              'strand': range_['strand']
                              })

    if check_uniq:
        if len(locations) > 1:
            raise RuntimeError('Sequence found in several subjects')

    return locations


if __name__ == '__main__':
    seq = 'AGACAAGTGGTGAAGAAKAAGATGATATGCAGCAATGCATTTCACCACTTTATATAGCATGGAGTGGATTTCTCCACCTCATTTAATAGTATGAAGTGGAGGCAGCCCCCCTCTACACCTGTCCACTAAGGCCAGCCCACAATCTGATCCCTTTTAATTTTTGCCTTGAGTGGTGGGGCCCATTGGATTAAATCAATCCAAATTAGCCAC'
    locations = locate_sequence(seq, config.TOMATO_GENOME_FASTA, check_uniq=False,
                                hsp_align_len_threshold=70, hsp_identity_threshold=95)
    pprint(locations)