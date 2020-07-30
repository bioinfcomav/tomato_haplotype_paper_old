
import config

from collections import Counter

import pandas

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import morphological
import passport
import pop_building
import plot
import labels


def compare_classifications(passports, morpho_classification, genetic_rank,
                            genetic_classes_to_ignore=None,
                            morpho_classes_to_ignore=None):

    if genetic_classes_to_ignore is None:
        genetic_classes_to_ignore = []
    if morpho_classes_to_ignore is None:
        morpho_classes_to_ignore = []

    genetic_classification = pop_building.get_classifications_for_classification_key_path(passports, genetic_rank)

    common_accs = sorted(set(genetic_classification.keys()).intersection(morpho_classification))

    common_morpho_classification = []
    common_genetic_classification = []
    final_accs = []
    for acc in common_accs:
        genetic_class = genetic_classification[acc]
        morpho_class = morpho_classification[acc]

        if genetic_class in genetic_classes_to_ignore:
            continue
        if morpho_class in morpho_classes_to_ignore:
            continue

        common_genetic_classification.append(labels.LABELS[genetic_class])
        common_morpho_classification.append(labels.LABELS[morpho_class])
        final_accs.append(acc)

    common_genetic_classification = pandas.Series(common_genetic_classification, index=final_accs)
    common_morpho_classification = pandas.Series(common_morpho_classification, index=final_accs)
    return {'genetic_classification': common_genetic_classification,
            'morpho_classification': common_morpho_classification}


if __name__ == '__main__':
    original_data = morphological.read_morphological_data()
    data = morphological.get_morphological_table_for_ordinal_traits()
    data.columns = [morphological.TRAIT_ABREVIATIONS[trait] for trait in data.columns]

    passports = passport.get_sample_passports()
    morpho_classification = morphological.read_morphological_classification()

    genetic_classes_to_ignore = [None, 'sp_x_sl', 'sp_x_sp']

    res = compare_classifications(passports, morpho_classification, config.RANK1,
                                  genetic_classes_to_ignore=genetic_classes_to_ignore,
                                  morpho_classes_to_ignore=[None, ''])

    plot_path = config.FIG_MORPHOLOGICAL_VS_MOLECULAR_CLASSIFICATION

    fig = Figure((10, 10))
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    #print(Counter(zip(res['genetic_classification'], res['morpho_classification'])))
    classes1_order = ['SP PE', 'SP Montane', 'SP PE Inter-Andean', 'SP EC', 'SLC EC', 'SLC CO',
                      'SLC MA', 'SLC PE', 'SLC world', 'SLL MX']
    classes2_order = ['SP PE', 'SP Montane', 'SP PE Inter-Andean', 'SP EC', 'SLC EC', 'SLC MA',
                      'SLC PE', 'SLL']

    plot.plot_table_classification_comparison(res['genetic_classification'],
                                              res['morpho_classification'],
                                              axes=axes,
                                              x_label='genetic',
                                              y_label='morphological',
                                              size_multiplier=10,
                                              classes1_order=classes1_order,
                                              classes2_order=classes2_order)


    fig.tight_layout()
    fig.savefig(plot_path)
