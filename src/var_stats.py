
import config

from variation.variations import VariationsH5

if False:
    vars_path = config.ORIG_H5
    variations = VariationsH5(str(vars_path), 'r')
    print('Orig vars')
    print('Num samples: ', len(variations.samples))
    # Num samples:  628
    print('Num vars: ', variations.num_variations)
    # Num vars:  25358552

if False:
    vars_path = config.TIER1_H5_LOWQ_085
    variations = VariationsH5(str(vars_path), 'r')
    print('tier1 lowq 85')
    print('Num samples: ', len(variations.samples))
    # Num samples:  598
    print('Num vars: ', variations.num_variations)
    # Num vars:  11803366

vars_path = config.WORKING_H5
variations = VariationsH5(str(vars_path), 'r')
print('working set')
# Num samples:  598
print('Num samples: ', len(variations.samples))
# Num vars:  33790
print('Num vars: ', variations.num_variations)
