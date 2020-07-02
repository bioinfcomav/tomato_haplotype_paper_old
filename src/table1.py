
import config

import csv

import bam_stats


def get_projects_for_samples(samples):
    projects = {}
    for sample in samples:
        if sample.startswith('TS-'):
            # https://doi.org/10.1016/j.cell.2017.12.019
            project = 'Zhu 2018'
        elif sample.startswith('BGV') or sample in ('Voyage', 'PAS014479', 'PI129033', 'PI129026', 'PI406890', 'CATIE-11106-1', 'PI129088', 'LA0767', 'LA2309', 'LA2697', 'PI487625', 'PI378994', 'Tegucigalpa', 'LA1712'):
            project = 'varitome'
        elif sample == 'heinz1706genome':
            project = 'The Tomato Genome Consortium 2012'
        elif sample.startswith('EA') or sample.startswith('TR'):
            # The 100 Tomato Genome Sequencing Consortium
            project = 'The 100 Tomato Genome Sequencing Consortium 2014'
        elif sample in ('Moneymaker', 'AlisaCraig'):
            project = 'The Tomato Genome Consortium 2012'
        elif sample in ('ECU1516B', 'ECU1032', 'ECU570B', 'GLP-12', 'GLP27B', 'GLP37A'):
            # tomato 6
            project = ''
        elif sample in ('LA0147', 'cervil', 'criollo', 'ferum', 'LA0147', 'LA1420', 'levovil', 'plovdiv', 'stupicke'):
            # DOI: 10.1186/1471-2164-14-791
            project = 'Causse 2013'
        elif sample == 'yellow_pear':
            # doi: 10.7717/peerj.793
            project = 'Strickler 2015'
        else:
            print(sample)
            assert False
        projects[sample] = project
    return projects


def write_mapping_stats(fhand):

    stats = bam_stats.read_bam_stats()['per_sample']

    projects = get_projects_for_samples(stats.keys())

    sorted_samples = sorted(stats.keys())
    sorted_samples = sorted(sorted_samples, key=lambda x: projects[x])

    writer = csv.writer(fhand)
    writer.writerow(['Accession', 'Reference', 'Mean cov.', 'Mean MAPQ57 cov.'])


    for sample in sorted_samples:
        sample_stats = stats[sample]
        writer.writerow([sample, projects[sample], sample_stats['coverages']['all_reads'], sample_stats['coverages']['mapq57']])


if __name__ == '__main__':
    path = config.FIGURES_DIR / 'mapping_stats.csv'
    write_mapping_stats(path.open('wt'))
