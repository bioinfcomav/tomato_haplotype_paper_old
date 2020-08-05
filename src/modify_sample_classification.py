samples_raw = '''
bgv008106	bgv008106	0.03127603	-0.0068651927	0.02580696
tr00019	tr00019	0.034147073	-0.015920116	0.0266611
bgv012640	bgv012640	0.09440451	0.08724301	-0.057934932
bgv008041	bgv008041	0.12135058	0.023373224	0.025726188
ts-66	ts-66	0.06790144	-0.07407922	0.13914981
bgv008189	bgv008189	0.109178886	-0.08627483	0.12605487
bgv007990	bgv007990	0.15870953	0.14256139	-0.049018674
ts-273	ts-273	0.184417	0.213392	-0.11925437
bgv007989	bgv007989	0.19446443	0.2114995	-0.11595837
bgv013945	bgv013945	0.07342253	0.07162432	0.047240395
ts-56	ts-56	0.17455785	0.13871354	-0.07824808
ts-242	ts-242	0.08856797	-0.006178598	0.067975
bgv008077	bgv008077	0.09444335	0.012995623	0.09520902
ts-74	ts-74	0.071566135	-0.03045697	0.11410013
ts-247	ts-247	0.088769265	0.023322979	0.07993729
bgv007981	bgv007981	0.106949024	-0.09192003	0.13310741
bgv008036	bgv008036	0.03282903	-0.0074944003	0.11588797
bgv008037	bgv008037	0.111176476	0.016282573	0.08093333

'''

rank = 'rank1'
classification = 'slc_pe_s'

lower_case_accession = True
ignore_missing_sample_error = False
add_missing_samples = True

import config
from classification import Classification

classifications = Classification(config.CLASSIFICATIONS, 'sample',
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
