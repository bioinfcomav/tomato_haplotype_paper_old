
import config

from collections import defaultdict


def get_classification_from_passport(classification_key_path, passport):

    passport_item = passport
    for key in classification_key_path:
        passport_item = passport_item.get(key, {})

    if not passport_item:
        passport_item = None

    if isinstance(passport_item, dict):
        print(passport_item)
        raise ValueError('The classification value should not be a dict but a str or None')

    classification = passport_item
    return classification


def get_classifications_for_classification_key_path(passports, classification_key_path,
                                                    classification_modifications=None):
    if classification_modifications is None:
        classification_modifications = {}

    classifications = {}
    for sample_id, passport in passports.items():
        classification = get_classification_from_passport(classification_key_path, passport)
        classification = classification_modifications.get(classification, classification)
        classifications[sample_id] = classification
    
    return classifications


def get_pops(pops_descriptions, passports, pops_to_merge=None):

    if pops_to_merge is None:
        pops_to_merge = {}

    provisional_classification_modifications = defaultdict(list)
    for merged_pop_name, original_pops in pops_to_merge.items():
        for original_pop in original_pops:
            provisional_classification_modifications[original_pop].append(merged_pop_name)

    classification_modifications = {}
    for original_pop, new_names in provisional_classification_modifications.items():
        if len(new_names) > 1:
            msg = f'{original_pop} is located within several pops to merge: '
            msg += ','.join(new_names)
            raise ValueError(msg)
        classification_modifications[original_pop] = new_names[0]

    classifications_by_key_path = {}
    for classification_key_path in pops_descriptions.keys():
        classifications = get_classifications_for_classification_key_path(passports,
                                                                          classification_key_path,
                                                                          classification_modifications)
        classifications_by_key_path[classification_key_path] = classifications

    pops = defaultdict(list)
    for classification_key_path, wanted_pops in pops_descriptions.items():
        classifications = classifications_by_key_path[classification_key_path]

        if wanted_pops == config.ALL_POPS:
            wanted_pops = set(classifications.values())

        for sample_id, pop in classifications.items():
            if pop not in wanted_pops:
                continue
            pops[pop].append(sample_id)

    pops = {pop: sorted(samples) for pop, samples in pops.items()}
    return pops


def get_merged_pops(pops_descriptions, passports):
    pass


if __name__ == '__main__':

    sample_passports = {'sample1': {'genetic_classification': {'rank1': 'pop1'}},
                        'sample2': {'genetic_classification': {'rank1': 'pop1'}},
                        'sample3': {'genetic_classification': {'rank1': 'pop2'}},
                        'sample4': {'genetic_classification': {'rank2': 'sub_pop1'}},
                        'sample5': {'morphological_classification': 'cultivated'}
                        }

    rank1 = ('genetic_classification', 'rank1', )
    rank2 = ('genetic_classification', 'rank2', )
    morphological = ('morphological_classification',)

    pops_descriptions = {rank1: ['pop1', 'pop2']}
    pops = get_pops(pops_descriptions, sample_passports)
    assert pops == {'pop1': ['sample1', 'sample2'], 'pop2': ['sample3']}

    pops_descriptions = {rank1: [None]}
    pops = get_pops(pops_descriptions, sample_passports)
    assert pops == {None: ['sample4', 'sample5']}

    pops_descriptions = {rank1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, sample_passports)
    assert pops == {'pop1': ['sample1', 'sample2'], 'pop2': ['sample3'], None: ['sample4', 'sample5']}

    pops_descriptions = {rank1: ['pop1'], rank2: ['sub_pop1']}
    pops = get_pops(pops_descriptions, sample_passports)
    assert pops == {'pop1': ['sample1', 'sample2'], 'sub_pop1': ['sample4']}

    pops_descriptions = {morphological: ['cultivated']}
    pops = get_pops(pops_descriptions, sample_passports)
    assert pops == {'cultivated': ['sample5']}

    pops_descriptions = {rank1: config.ALL_POPS}
    try:
        pops = get_pops(pops_descriptions, sample_passports,
                        pops_to_merge={'merged_pop':['pop1', 'pop2'], 'pop1': ['pop1']})
        assert RuntimeError('ValueError expected here')
    except ValueError:
        pass

    pops_descriptions = {rank1: config.ALL_POPS}
    pops = get_pops(pops_descriptions, sample_passports,
                    pops_to_merge={'merged_pop':['pop1', 'pop2'] }  )
    assert pops == {'merged_pop': ['sample1', 'sample2', 'sample3'], None: ['sample4', 'sample5']}
