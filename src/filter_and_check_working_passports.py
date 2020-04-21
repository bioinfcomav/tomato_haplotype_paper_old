import config

import csv
from pandas import read_csv, isnull

from variation.variations import VariationsH5


if __name__ == "__main__":
    passports_csv = read_csv(config.PASSPORTS_CSV)
    passports = {}
    for passport in passports_csv.iterrows():
        passports[passport[1]['accession']] = passport[1]

    variations = VariationsH5(str(config.WORKING_H5), 'r')
    samples = variations.samples

    fieldnames = ['accession', 'country', 'latitude', 'longitude', 'species', 'variety', 'region', 'BGV', 'LA', 'EA', 'PI', 'LYC', 'CATORT', 'CGN']
    writer = csv.DictWriter(str(config.WORKING_PASSPORTS_CSV), delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    for sample in sorted(samples):
        row = {}
        for field in fieldnames:
            value = passports[sample][field]
            if isnull(value):
                value = ''
            row[field] = value
        writer.writerow(row)
