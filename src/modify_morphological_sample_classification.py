samples_raw = '''

bgv006231

'''

rank = 'morpho_class'
classification = 'slc_ecu'

lower_case_accession = True
ignore_missing_sample_error = False
add_missing_samples = True

import config
from classification import Classification

classifications = Classification(config.CLASSIFICATIONS_MORPHOLOGICAL, 'sample',
                                history_dir=config.CLASSIFICATION_HISTORY_DIR,
                                ignore_missing_sample_error=ignore_missing_sample_error,
                                add_missing_samples=add_missing_samples)

samples_to_change = []

for line in samples_raw.splitlines():
    line = line.strip()
    if not line:
        continue
    if line.startswith('Colour'):
        continue
    accession = line.split('\t')[0]
    if lower_case_accession:
        accession = accession.lower()
    samples_to_change.append(accession)

print(samples_to_change)
print(rank, classification)
classifications.modify_samples_classification(samples_to_change,
                                              rank,
                                              classification)
classifications.save()
