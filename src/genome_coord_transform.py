
from collections import OrderedDict

TOMATO_GENOME_FAI = '''SL4.0ch01       90863682        9803994 60      61
SL4.0ch02       53473368        102182082       60      61
SL4.0ch03       65298490        156546685       60      61
SL4.0ch04       64459972        222933495       60      61
SL4.0ch05       65269487        288467812       60      61
SL4.0ch06       47258699        354825136       60      61
SL4.0ch07       67883646        402871491       60      61
SL4.0ch08       63995357        471886544       60      61
SL4.0ch09       68513564        536948503       60      61
SL4.0ch10       64792705        606603971       60      61
SL4.0ch11       54379777        672476567       60      61
SL4.0ch12       66688036        727762686       60      61
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


if __name__ == '__main__':
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
