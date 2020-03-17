
import hashlib
import pickle

from variation.variations.distance import calc_pairwise_distance


def calc_kosman_dists(variations, method='kosman', cache_dir=None):
    if cache_dir:
        key = ','.join(sorted(variations.samples)) + str(variations.num_variations) + method
        key = hashlib.md5(key.encode()).hexdigest()
        cache_path = cache_dir / ('dists' + key + 'dists')
        if cache_path.exists():
            return pickle.load(cache_path.open('rb'))

    print('Calculating distance')
    dists = calc_pairwise_distance(variations, chunk_size=1000, method=method)
    result = {'dists': dists, 'samples': variations.samples[:]}

    if cache_dir:
        pickle.dump(result, cache_path.open('wb'))

    return result
