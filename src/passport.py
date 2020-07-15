
from csv import DictReader
import config

from pprint import pprint
import csv


def get_acc_passports():
    reader = csv.DictReader(config.PASSPORTS_CSV.open('rt'), delimiter=',')

    accessions = {}
    for row in reader:

        row = {key.strip(): value.strip() for key, value in row.items()}

        accession_id = row['accession']
        accession_info = {}
        if row['country']:
            accession_info['country'] = row['country']
        if row['species']:
            taxon = row['species']
            if row['variety']:
                taxon += f' var. {row["variety"]}'
            accession_info['taxon'] = taxon

        if row['latitude']:
            accession_info['latitude'] = float(row['latitude'])
            accession_info['longitude'] = float(row['longitude'])

        if row['region']:
            accession_info['geo_region'] = row['region']

        accessions[accession_id] = accession_info
    return accessions


def get_classifications():
    reader = csv.DictReader(config.CLASSIFICATIONS.open('rt'))
    classifications = {}

    for row in reader:
        classification = {rank: row[rank] if row[rank] else None for rank in config.CLASSIFICATION_RANKS}
        classifications[row['sample']] = classification
    return classifications


def get_razifard_classifications():

    reader = DictReader(config.RAZIFARD_CLASSIFICATION.open('rt'), delimiter='\t')

    classification = {}
    for sample in reader:
        if sample['Razifard_pop']:
            classification[sample['sample']] = sample['Razifard_pop']
    return classification


def get_sample_passports():
    
    accessions = get_acc_passports()

    # for this project there is only one sample per accession
    samples = accessions

    classifications = get_classifications()

    razifard_classifications = get_razifard_classifications()
    
    for sample, sample_info in samples.items():
        classification = classifications.get(sample)
        sample_info['classification'] = classification
        if razifard_classifications.get(sample):
            sample_info['razifard_classification'] = razifard_classifications.get(sample)
    return samples


if __name__ == '__main__':
    get_sample_passports()
