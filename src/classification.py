import csv
from pprint import pprint
from collections import OrderedDict
import time
import shutil
import pathlib

import config


class Classification:
    def __init__(self, csv_path, sample_col_name, delimiter=',', history_dir=None,
                ignore_missing_sample_error=False, add_missing_samples=False):
        self.csv_path = csv_path
        self.sample_col_name = sample_col_name
        self.delimiter = delimiter
        self._classification = None
        self._ranks = None
        self.history_dir = history_dir
        self.ignore_missing_sample_error = ignore_missing_sample_error
        self.add_missing_samples = add_missing_samples
        self._load()

    def _load(self):
        reader = csv.DictReader(self.csv_path.open('rt'),
                                delimiter=self.delimiter)
        self._ranks = [rank for rank in reader.fieldnames if rank != self.sample_col_name]
        self._classification = OrderedDict()
        for row in reader:
            sample_classification = {rank: row[rank] for rank in self._ranks}
            self._classification[row[self.sample_col_name]] = sample_classification

    def save(self):
        path = self.csv_path
        fhand = self.csv_path.open('wt')
        writer = csv.DictWriter(fhand,
                                fieldnames=[self.sample_col_name] + self._ranks)
        writer.writeheader()
        for sample_id, sample_classification in self._classification.items():
            row = sample_classification.copy()
            row[self.sample_col_name] = sample_id
            writer.writerow(row)
        fhand.flush()

        if self.history_dir:
            datestamp = time.strftime("%Y-%m-%d_%H:%M:%S", time.gmtime())
            history_path = self.history_dir / f'{path.stem}_{datestamp}.csv'
            shutil.copyfile(fhand.name, history_path)

    def modify_sample_classification(self, sample_id, rank, new_classification):
        try:
            self._classification[sample_id][rank] = new_classification
        except KeyError:
            if self.ignore_missing_sample_error:
                pass
            elif self.add_missing_samples:
                self._classification[sample_id] = {rank: new_classification}
            else:
                raise

    def modify_samples_classification(self, sample_ids, rank, new_classification,
                                      ignore_missing_sample_error=False):
        for sample_id in sample_ids:
            self.modify_sample_classification(sample_id, rank, new_classification)

    def update_classification_names(self, rank, old_name, new_name):
        for sample_id, sample_classification in self._classification.items():
            if sample_classification[rank] == old_name:
                self.modify_sample_classification(sample_id, rank, new_name)

if __name__ == '__main__':
    classification = Classification(config.CLASSIFICATIONS, 'sample',
                                    history_dir=config.CLASSIFICATION_HISTORY_DIR)
    classification.save()
