
import config

import pathlib
from pprint import pprint
from collections import defaultdict

import numpy

import genome_coord_transform


def read_bam_list():
    sample_names = [pathlib.Path(line).stem for line in config.BAM_LIST.open('rt')]
    return {'bam_stem_names': sample_names}
    

def _read_bam_stats(path):

    coverages = {}
    num_mapped_bases = None
    for line in path.open('rt'):
        #print(line)
        if line.startswith('COV'):
            items = line.split()
            coverages[int(items[2])] = int(items[3])
        if line.startswith('SN'):
            if line.startswith('SN	bases mapped (cigar):'):
                num_mapped_bases = int(line.split()[-4])
            if line.startswith('SN	average quality:'):
                average_quality = float(line.split()[-1])

    res = {'coverages': coverages,
           'average_quality': average_quality,
           'num_mapped_bases': num_mapped_bases}

    return res


def _get_bamstat_path_mapq57(sample):
    bamstats_path1 = config.BAM_STATS_DIR / (sample + '.bamstats.mapq57.txt')
    if bamstats_path1.exists():
        return bamstats_path1
    raise RuntimeError(f'bamstat file not found: {bamstats_path1}')


def _get_header_path_for_bam_stem(bam_stem_name):
    path = config.BAM_STATS_DIR / (bam_stem_name + '.bam.header')

    if not path.exists():
        raise RuntimeError(f'header file missing: {path}')
    return path


def _get_sample_for_bam(bam_stem_name):
    path = _get_header_path_for_bam_stem(bam_stem_name)

    for line in path.open('rt'):
        if line.startswith('@RG'):
            rg = dict([item.split(':') for item in line[3:].split()])
            sample = rg['SM']
    return sample


def _get_genome_size(bam_stem_name):
    path = _get_header_path_for_bam_stem(bam_stem_name)

    size = 0
    for line in path.open('rt'):
        if line.startswith('@SQ'):
            sq = dict([item.split(':') for item in line[3:].split()])
            size += int(sq['LN'])
    return size


def read_bam_stats():
    bam_stem_names = read_bam_list()['bam_stem_names']

    results = {}
    genome_sizes = defaultdict(dict)
    for bam_stem_name in bam_stem_names:

        sample = _get_sample_for_bam(bam_stem_name)
        if sample not in results:
            results[sample] = {}
        results[sample][bam_stem_name] = {}
    
        bamstats_path = config.BAM_STATS_DIR / (bam_stem_name + '.bamstats.txt')
        try:
            results[sample][bam_stem_name]['all_reads'] = _read_bam_stats(bamstats_path)
        except (FileNotFoundError, RuntimeError):
            raise

        genome_sizes[sample][bam_stem_name] = _get_genome_size(bam_stem_name)

        bamstats_path = _get_bamstat_path_mapq57(bam_stem_name)
        try:
            results[sample][bam_stem_name]['mapq57'] = _read_bam_stats(bamstats_path)
        except (FileNotFoundError, RuntimeError):
            raise

    kind_map_filterings = {kind_map_filtering for sample_stats in results.values() for bam_stats in sample_stats.values() for kind_map_filtering in bam_stats.keys()}

    results_per_sample = {}
    for sample, sample_stats in results.items():
        num_bases_mapped = {kind_map_filtering: 0 for kind_map_filtering in kind_map_filterings}
        genome_sizes_for_sample = []
        for bam_stem_name, bam_stats in sample_stats.items():
            genome_sizes_for_sample.append(genome_sizes[sample][bam_stem_name])
            for kind_map_filtering in kind_map_filterings:
                num_bases_mapped[kind_map_filtering] += bam_stats.get(kind_map_filtering, {}).get('num_mapped_bases', 0)
    
        assert all([genome_sizes_for_sample[0] == size for size in genome_sizes_for_sample])
        genome_size = genome_sizes_for_sample[0]

        results_per_sample[sample] = {'num_bases_mapped': num_bases_mapped}

        coverages = {kind_map_filtering: num_bases_mapped / genome_size for kind_map_filtering, num_bases_mapped in results_per_sample[sample]['num_bases_mapped'].items()}
        results_per_sample[sample]['coverages'] = coverages

    return {'per_bam': results, 'per_sample': results_per_sample}


if __name__ == '__main__':
    res = read_bam_stats()

