samples_raw = '''
bgv006825	bgv006825	2.3738673	0.6126412	-1.6645778
bgv013945	bgv013945	2.5002782	1.2262374	-1.4652416
bgv008077	bgv008077	3.1624382	1.3311117	-1.8544629
bgv012625	bgv012625	2.6388392	0.63737863	-1.5750884
bgv007981	bgv007981	3.6100519	1.1890577	-2.190068
bgv007936	bgv007936	2.9165926	1.4236209	-1.4269761
bgv007870	bgv007870	3.3410258	1.4118094	-1.1278807
bgv007873	bgv007873	3.733647	1.7947648	-1.5208321
bgv007865	bgv007865	3.584294	1.8609052	-1.667412

'''

rank = 'morpho_class'
classification = 'sll'

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
