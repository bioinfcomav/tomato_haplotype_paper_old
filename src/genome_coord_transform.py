
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


def get_genome_sizes():
    lengths = {}
    for line in TOMATO_GENOME_FAI.splitlines():
        if not line.strip():
            continue
        items = line.split()
        chrom = items[0]
        length = int(items[1])
        lengths[chrom] = length
    return lengths


class GenomeCoordinateConverter:
    def __init__(self, chrom_lens=None):
        if chrom_lens is None:
            chrom_lens = get_genome_sizes()

        sorted_chroms = sorted(chrom_lens.keys())
        offsets = OrderedDict()
        offset = 0
        chrom_spans = {}
        for chrom in sorted_chroms:
            chrom_len = chrom_lens[chrom]
            chrom_name = chrom.lower()
            offsets[chrom_name] = offset
            chrom_start = offset
            offset += chrom_len
            chrom_end = offset
            self.genome_size = offset
            chrom_spans[chrom_name] = (chrom_start, chrom_end)
        self._offsets = offsets
        self.chrom_spans = chrom_spans
        self.chrom_lens = chrom_lens

    def transform_coordinate(self, chrom, pos):
        if isinstance(chrom, bytes):
            chrom = chrom.decode()
        try:
            return self._offsets[chrom.lower()] + pos
        except KeyError:
            raise ValueError('Unknown chromosome: ' + chrom)


class PositionInPericentromericRegion(Exception):
    pass


class GenomeCoordinateConverter2:
    def __init__(self):
        res = get_euchromatic_regions()
        self._euchromatic_regions = res['euchromatic_regions']
        self.chrom_lens = res['euchromatic_sizes']

        self.sorted_chroms = sorted(self._euchromatic_regions.keys())
        euchromatic_regions = self._euchromatic_regions

        chrom_offsets = {}
        chrom_spans = {}
        pericentromeric_starts = {}
        offset = 0
        for chrom in self.sorted_chroms:
            chrom_offsets[chrom] = offset

            chrom_start = offset
            offset += euchromatic_regions[chrom][0][1] + euchromatic_regions[chrom][1][1] - euchromatic_regions[chrom][1][0]
            pericentromeric_starts[chrom] = euchromatic_regions[chrom][0][1]
            chrom_end = offset
            self.genome_size = offset
            chrom_spans[chrom] = (chrom_start, chrom_end)

        self.chrom_offsets = chrom_offsets
        self.chrom_spans = chrom_spans
        self.pericentromeric_starts = pericentromeric_starts

    def _get_offset(self, chrom, pos):
        try:
            chrom_offset = self.chrom_offsets[chrom]
        except KeyError:
            raise KeyError('Unknown chromosome')

        euchromatic_regions = self._euchromatic_regions[chrom]
        if pos <= euchromatic_regions[0][1]:
            return chrom_offset, 0
        if euchromatic_regions[0][1] < pos and pos < euchromatic_regions[1][0]:
            raise PositionInPericentromericRegion()
        if pos >= euchromatic_regions[1][0]:
            return chrom_offset + euchromatic_regions[0][1], euchromatic_regions[1][0]

    def transform_coordinate(self, chrom, pos):
        if isinstance(chrom, bytes):
            chrom = chrom.decode()
        offset, to_remove_from_pos = self._get_offset(chrom, pos)
        return offset + pos - to_remove_from_pos


def get_euchromatic_regions():
    pericentromeric_regions = {}
    for line in PERICENTROMERIC_REGIONS_BED.splitlines():
        chrom, peri_start, peri_end = line.split()
        peri_start = int(peri_start)
        if peri_start < 1:
             peri_start = 2
        peri_end = int(peri_end)
        pericentromeric_regions[chrom] = (peri_start, peri_end)

    chrom_sizes = get_genome_sizes()
    euchromatin_regions = {}
    euchromatic_sizes = {}
    tuples = []
    for chrom, size in chrom_sizes.items():
        peri_region = pericentromeric_regions[chrom]
        assert peri_region[1] < size
        if peri_region[0] > 10:
            tuples.append((chrom, 1, peri_region[0]))
        tuples.append((chrom, peri_region[1], size))
        euchromatin_regions[chrom] = (1, peri_region[0]), (peri_region[1], size)
        euchromatic_sizes[chrom] = (peri_region[0] - 1) + (size - peri_region[1])

    return {'euchromatic_regions': euchromatin_regions,
            'euchromatic_sizes': euchromatic_sizes,
            'euchromatic_regions_as_tuples': tuples}



if __name__ == '__main__':
    get_euchromatic_regions()
    coord_converter = GenomeCoordinateConverter2()
    print(coord_converter.transform_coordinate('SL2.50ch01', 10))
    try:
        coord_converter.transform_coordinate('SL2.50ch01', 5488554)
    except PositionInPericentromericRegion:
        pass
    print(coord_converter.transform_coordinate('SL2.50ch01', 74024605))

    coord_converter = GenomeCoordinateConverter()
    print(coord_converter.transform_coordinate('SL2.50ch01', 98543443))
    print(coord_converter.transform_coordinate('SL2.50ch02', 2))
