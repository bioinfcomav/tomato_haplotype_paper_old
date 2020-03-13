
from collections import defaultdict

KEEP = 'keep'
REMOVE = 'remove'


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
            keys, pops, action = criterion
            passport_item = passport
            for key in keys:
                passport_item = passport_item[key]
            pop = passport_item
            pops_seen[keys].add(pop)

            if pop in pops:
                if action == KEEP:
                    samples_to_keep.add(sample)
                if action == REMOVE:
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