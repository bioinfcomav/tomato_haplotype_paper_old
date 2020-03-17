
from collections import OrderedDict

TOMATO_GENOME_FAI = '''SL2.50ch01	98543444	22117358	70	71
SL2.50ch02	55340444	122068579	70	71
SL2.50ch03	70787664	178199614	70	71
SL2.50ch04	66470942	249998544	70	71
SL2.50ch05	65875088	317419084	70	71
SL2.50ch06	49751636	384235258	70	71
SL2.50ch07	68045021	434697645	70	71
SL2.50ch08	65866657	503714751	70	71
SL2.50ch09	72482091	570522374	70	71
SL2.50ch10	65527505	644039937	70	71
SL2.50ch11	56302525	710503563	70	71
SL2.50ch12	67145203	767610423	70	71
'''

PERICENTROMERIC_REGIONS_BED = '''SL2.50ch01	5488553	74024603
SL2.50ch02	0	30493730
SL2.50ch03	16493431	50407653
SL2.50ch04	7406888	50551374
SL2.50ch05	9881466	58473554
SL2.50ch06	3861081	33077717
SL2.50ch07	4056987	58629226
SL2.50ch08	4670213	54625578
SL2.50ch09	6225214	63773642
SL2.50ch10	3775719	55840828
SL2.50ch11	10947270	48379978
SL2.50ch12	5879033	61255621
'''


def get_chrom_sizes():
    lengths = OrderedDict()
    for line in TOMATO_GENOME_FAI.splitlines():
        if not line.strip():
            continue
        items = line.split()
        chrom = items[0]
        length = int(items[1])
        lengths[chrom] = length
    return lengths


def get_pericentromic_regions():
    regions = OrderedDict()
    for line in PERICENTROMERIC_REGIONS_BED.splitlines():
        if not line.strip():
            continue
        items = line.split()
        chrom = items[0]
        start = int(items[1])
        end = int(items[2])
        regions[chrom] = (start, end)
    return regions


if __name__ == '__main__':
    print(get_chrom_sizes())
    print(get_pericentromic_regions())