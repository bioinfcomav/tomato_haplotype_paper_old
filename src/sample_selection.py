import config

from collections import defaultdict

from pop_building import get_classification_from_passport


def get_samples_for_criteria(all_samples, sample_passports, criteria,
                             skip_samples_with_no_passport=False):
    samples_to_keep = criteria.get('samples_to_keep', [])
    samples_to_remove = criteria.get('samples_to_remove', [])
    samples_to_keep = set(samples_to_keep)
    samples_to_remove = set(samples_to_remove)

    criteria = criteria.get('criteria', [])

    if not criteria and not samples_to_keep and not samples_to_remove:
        return all_samples

    samples_to_use = []
    pops_seen = defaultdict(set)
    for sample in all_samples:
        if sample in samples_to_keep:
            samples_to_use.append(sample)
        if sample in samples_to_remove:
            continue

        try:
            passport = sample_passports[sample]
        except KeyError:
            if skip_samples_with_no_passport:
                continue
            raise

        keep_it = False
        remove_it = False
        for criterion in criteria:
            classification_key_path, pops, action = criterion
            pop = get_classification_from_passport(classification_key_path, passport)

            pops_seen[classification_key_path].add(pop)

            if pop in pops:
                if action == config.KEEP:
                    samples_to_keep.add(sample)
                if action == config.REMOVE:
                    samples_to_remove.add(sample)

    if samples_to_keep:
        samples = [sample for sample in all_samples if sample in samples_to_keep]
    else:
        samples = all_samples

    samples = [sample for sample in samples if sample not in samples_to_remove]

    if not samples:
        print('Populations avaliable:')
        for classification_keys, pops in pops_seen.items():
            print(classification_keys)
            for pop in sorted(pops, key=lambda x: str(x)):
                print(f'\t{pop}')

        raise ValueError('No samples left for the given criteria')

    return samples

if __name__ == '__main__':

    passports = {'sample1': {'genetic_classification': {'rank1': 'pop1'}},
                 'sample2': {'genetic_classification': {'rank1': 'pop1'}},
                 'sample3': {'genetic_classification': {'rank1': 'pop2'}},
                 'sample4': {'genetic_classification': {'rank2': 'sub_pop1'}},
                 'sample5': {'morphological_classification': 'cultivated'}
                        }

    rank1 = 'genetic_classification', 'rank1'
    criteria = {'criteria': [(rank1, ['pop1'], KEEP)]}

    all_samples = list(passports.keys())

    samples_to_use = get_samples_for_criteria(all_samples,
                                              passports,
                                              criteria,
                                              skip_samples_with_no_passport=True)

    assert sorted(samples_to_use) == ['sample1', 'sample2']