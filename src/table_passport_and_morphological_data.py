
import config

import morphological

if __name__ == '__main__':
    original_data = morphological.read_morphological_data()

    morpho_classification = morphological.read_morphological_classification()

    path = config.TABLE_PASSPORT_AND_MORPHOLOGICAL_DATA

    morphological.write_morpho_csv(path, original_data, morpho_classification, write_collection_id=False)
