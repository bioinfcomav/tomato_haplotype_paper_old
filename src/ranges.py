
from pprint import pprint

import pandas

import pyranges


def merge_ranges(ranges):
    ranges = create_pyranges_df(ranges)

    use_strand = 'Strand' in ranges.df.columns

    ranges = ranges.merge(count=False, strand=use_strand)
    return pyranges_to_dict_ranges(ranges)


def pyranges_to_dict_ranges(ranges):
    dict_ranges = []

    use_strand = 'Strand' in ranges.df.columns

    if use_strand:
        for chrom, start, end, strand in zip(ranges.Chromosome, ranges.Start, ranges.End, ranges.Strand):
            dict_ranges.append({'chrom': chrom, 'start': start, 'end': end, 'strand': strand})
    else:
        for chrom, start, end in zip(ranges.Chromosome, ranges.Start, ranges.End):
            dict_ranges.append({'chrom': chrom, 'start': start, 'end': end})

    return dict_ranges


def create_pyranges_df(ranges):
    chroms = []
    starts = []
    ends = []
    strands = []
    for range_ in ranges:
        chroms.append(range_['chrom'])
        starts.append(range_['start'])
        ends.append(range_['end'])
        try:
            strands.append(range_['strand'])
        except:
            continue
    if strands and len(strands) != len(chroms):
        raise ValueError('There are strands only for some ranges')

    dict_ = {'Chromosome': chroms,
             'Start': starts,
             'End': ends}
    if strands:
        dict_['Strand'] = strands

    ranges_df = pandas.DataFrame(dict_)
    ranges = pyranges.PyRanges(ranges_df)
    return ranges
    

if __name__ == '__main__':
    ranges = [{'chrom': '1', 'start': 1, 'end': 10}, {'chrom': '1', 'start': 5, 'end': 15}]
    merge_ranges(ranges)

    ranges = [{'chrom': '1', 'start': 1, 'end': 10, 'strand': '+'}, {'chrom': '1', 'start': 5, 'end': 15, 'strand': '-'}]
    print(merge_ranges(ranges))
