samples_raw = '''
bgv008106	bgv008106	2.67129	1.1873915	-0.31173804
bgv008100	bgv008100	3.6040657	0.31009722	0.73302287
bgv006852	bgv006852	3.1033254	0.88562113	-0.31927118
bgv013946	bgv013946	4.51641	1.2837021	0.07271313
pi-129088	pi-129088	2.8838112	1.2237607	-0.53289783
mex-059	mex-059	2.8927817	0.5653325	-0.38451248
bgv008354	bgv008354	2.6646264	0.13605328	0.00489308
bgv008095	bgv008095	4.10329	1.1181567	-0.040323567
bgv008098	bgv008098	3.4253256	1.0420372	-0.056709263
bgv006931	bgv006931	3.0383008	1.1615485	-0.09333816
bgv006934	bgv006934	3.2966845	0.49882564	0.02679204
bgv016056	bgv016056	3.0074892	-0.2641853	0.4598421
bgv016055	bgv016055	3.400523	0.1810674	0.07334475
bgv016054	bgv016054	2.803849	0.45010203	-0.097947784
bgv016053	bgv016053	4.1694856	1.0298377	0.35384727
bgv016052	bgv016052	2.8569617	0.32244298	-0.0030153382
bgv016050	bgv016050	2.9060152	0.51368505	0.21710311
bgv005912	bgv005912	2.7916234	0.47552207	0.102078564
bgv008060	bgv008060	3.6131747	0.35659155	1.7055478
bgv008065	bgv008065	3.7214243	-0.036163967	0.39710784
bgv015729	bgv015729	2.8338711	1.1243178	0.5014068
bgv015726	bgv015726	3.7876198	-0.124482855	0.79127866
bgv015727	bgv015727	3.0233738	0.16742386	-0.11881188
bgv015725	bgv015725	3.04462	0.019863617	1.2353168
bgv013156	bgv013156	3.050049	0.46942878	-0.46845543
bgv004584	bgv004584	2.8838112	1.2237607	-0.53289783
bgv007982	bgv007982	4.253848	0.16606647	0.8467838

'''

rank = 'morpho_class'
classification = 'slc_big'

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
