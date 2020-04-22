import config

import csv
from pandas import read_csv, isnull

from variation.variations import VariationsH5


if __name__ == "__main__":
    passports_csv = read_csv(config.CLASSIFICATIONS)
    passports = {}
    for passport in passports_csv.iterrows():
        passports[passport[1]['sample']] = passport[1]

    variations = VariationsH5(str(config.WORKING_H5), 'r')
    samples = variations.samples

    fieldnames = ['sample', 'rank1', 'rank2', 'morpho_type']
    writer = csv.DictWriter(open(str(config.WORKING_CLASSIFICATION_CSV), 'w'), fieldnames=fieldnames)
    writer.writeheader()

    for sample in sorted(samples):
        row = {}
        for field in fieldnames:
            value = passports[sample][field]
            if isnull(value):
                value = ''
            row[field] = value
        writer.writerow(row)
