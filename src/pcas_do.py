
import config

from pprint import pprint

from variation.variations import VariationsH5

from passport import get_sample_passports
from sample_selection import get_samples_for_criteria, KEEP, REMOVE
from snp_filtering import filter_variations
import colors


RANK1 = 'classification', 'rank1'
RANK2 = 'classification', 'rank2'



def get_sample_selection_criteria(all_samples):

    criteria = []
    samples_to_remove = []
    samples_to_keep = ['hola']

    if True:
        criteria.append((RANK1, ['sp_pe', 'sp_ec', 'sp_inter-andean'], KEEP))
        #criteria.append((RANK2, ['slc_ecu_big'], REMOVE))

    return {'criteria': criteria, 'samples_to_remove': samples_to_remove,
            'samples_to_keep': samples_to_keep}


if __name__ == '__main__':

    debug = True

    filter_by_maf = True
    chunk_size = 1000

    if debug:
        max_chunks_to_process = None
        cache_dir=config.CACHE_DIR
    else:
        max_chunks_to_process = 10
        cache_dir = None


    passports = get_sample_passports()

    vars_path = config.TIER1_PHASED_LOW_QUAL_09_MISSING_085
    variations = VariationsH5(str(vars_path), 'r')

    all_samples = variations.samples
    print(variations.num_variations)

    criteria = get_sample_selection_criteria(all_samples)

    samples_to_use = get_samples_for_criteria(all_samples,
                                              passports,
                                              criteria,
                                              skip_samples_with_no_passport=True)
    
    if filter_by_maf:
        print(colors.TERMINAL_BLUE + 'Doing a MAF filtering you are reducing the distance of the samples that do not belong to a well represented population' + colors.TERMINAL_ENDC)
        max_maf = config.TIER2_PCA_MAX_MAF
    else:
        max_maf = None

    print(variations.num_variations)
    tier2_vars = filter_variations(variations,
                                   chunk_size=chunk_size,
                                   samples_to_keep=samples_to_use,
                                   cache_dir=cache_dir,
                                   max_mac=config.TIER2_MAX_MAC,
                                   max_maf=max_maf,
                                   min_called=config.TIER2_MIN_CALLED,
                                   max_het=config.TIER2_MAX_HET,
                                   min_call_for_het=config.TIER2_MAX_HET_MIN_CALL_DP,
                                   kept_fields=config.RELEVANT_FIELDS,
                                   max_chunks_to_process=max_chunks_to_process,
                                   remove_non_variable_snvs=True,
                                   verbose=True
                                   )

    print(tier2_vars.num_variations)