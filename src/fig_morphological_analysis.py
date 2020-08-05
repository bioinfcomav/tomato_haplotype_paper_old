
import config

import pandas

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import morphological
import passport
import labels
import colors
import matplotlib_support
import plot
import morphological_progression
import pop_building


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
        common_morpho_classification.append(labels.LABELS.get(morpho_class, morpho_class))
        final_accs.append(acc)

    common_genetic_classification = pandas.Series(common_genetic_classification, index=final_accs)
    common_morpho_classification = pandas.Series(common_morpho_classification, index=final_accs)
    return {'genetic_classification': common_genetic_classification,
            'morpho_classification': common_morpho_classification}


if __name__ == '__main__':
    original_data = morphological.read_morphological_data()
    data = morphological.get_morphological_table_for_ordinal_traits()
    data.columns = [morphological.TRAIT_ABREVIATIONS[trait] for trait in data.columns]
    #print(Counter([acc_data['characterization']['stripy_fruit'] for acc_data in read_morphological_data().values()]))

    data_no_na = morphological.fill_missing_data_with_means(data, max_missing_values=8)
    pca_result = morphological.do_pca(data_no_na, num_dims=3)

    passports = passport.get_sample_passports()
    morpho_classification = morphological.read_morphological_classification()

    morpho_classification = {acc: labels.LABELS[klass] for acc, klass in morpho_classification.items()}

    plot_path = config.FIG_MORPHOLOGICAL_ANALYSIS

    fig = Figure((10, 10))
    FigureCanvas(fig) # Don't remove it or savefig will fail later

    axes = matplotlib_support.add_axes(fig, row_idx=0, col_idx=0,
                                       axes_col_widths=[0.5, .5], axes_row_heights=[0.5, 0.5],
                                       bottom_margin=0.1, left_margin=0.13)
    color_schema = colors.ColorSchema(colors.CLASSIFICATION_RANK1_COLORS)
    morphological.plot_morphological_pca(pca_result, axes, morpho_classification,
                                         color_schema=color_schema,
                                         plot_principal_components=True,
                                         principal_component_scaler=3)
    matplotlib_support.set_axes_background(axes)

    axes = matplotlib_support.add_axes(fig, row_idx=0, col_idx=1,
                                       axes_col_widths=[0.5, .5], axes_row_heights=[0.5, 0.5],
                                       left_margin=0.20, bottom_margin=0.25)

    genetic_classes_to_ignore = [None, 'sp_x_sl', 'slc_co']
    morpho_classes_to_ignore = [None, '', 'Unclassified']
    #genetic_classes_to_ignore = [None]

    res = compare_classifications(passports, morpho_classification, config.RANK1,
                                  genetic_classes_to_ignore=genetic_classes_to_ignore,
                                  morpho_classes_to_ignore=morpho_classes_to_ignore)

    # genetic
    classes1_order = ['SP Pe', 'SP Montane', 'SP Ec', 'SP x SP', 'SLC Ec', 'SLC Co',
                      'SLC MA', 'SLC Pe N', 'SLC Pe S', 'SLC world', 'SLL MX', 'SP x SL', 'SP-SL', 'Unclassified']

    classes2_order = ['SP Pe', 'SP Montane', 'SP Ec', 'SLC Ec',
                      'SLC small', 'SLC big', 'SLL', 'SP x SL', 'SP-SL', 'Unclassified']

    plot.plot_table_classification_comparison(res['genetic_classification'],
                                              res['morpho_classification'],
                                              axes=axes,
                                              x_label='genetic',
                                              y_label='morphological',
                                              size_multiplier=10,
                                              classes1_order=classes1_order,
                                              classes2_order=classes2_order)

    morpho_classes = ['sp_pe', 'sp_montane', 'sp_ec',
                      'slc_ec', 'slc_small', 'slc_big', 'sll']
    morpho_classes = [labels.LABELS[klass] for klass in morpho_classes]
    background_colors = {klass: colors.modify_color(color_schema[klass], saturation_mod=-0.1, luminosity_mod=0.15) for klass in morpho_classes}

    inflorescence_traits = ['presence_of_irregular_inflorescence', 'style_exsertion', 'style_curvature',
              'petal_position', 'petal_width', 'inflorescence_ordinal']
    stem_traits = ['stem_width', 'stem_hairiness']
    leaf_traits = ['leaf_type', 'leaflet_margin']
    fruit_traits = ['ribbing', 'fruit_size', 'peanut_fruit_shape', 'stripy_fruit']

    all_traits = (('Inflorescence', inflorescence_traits),
                  ('Stem', stem_traits),
                  ('Leaf', leaf_traits),
                  ('Fruit', fruit_traits) 
                 )
    for idx, (structure, traits) in enumerate(all_traits):
        traits = [morphological.TRAIT_ABREVIATIONS[trait] for trait in traits]
        row_idx = idx // 2 + 1
        col_idx = idx % 2
        top_margin = 0
        left_margin = 0.12 if col_idx == 0 else 0.05
        bottom_margin = 0.05 if row_idx == 1 else 0.25
        axes = matplotlib_support.add_axes(fig, row_idx=row_idx, col_idx=col_idx,
                                           top_margin=top_margin, left_margin=left_margin, bottom_margin=bottom_margin,
                                           axes_col_widths=[0.5, .5], axes_row_heights=[0.5, 0.22, 0.29])
        morphological_progression.plot_morphological_progression(data, morpho_classification, axes, sorted_morpho_classes=morpho_classes,
                                                                traits=traits,
                                                                background_colors=background_colors,
                                                                bar_alpha=0.4)
        if row_idx == 1:
            matplotlib_support.turn_off_x_axis(axes)
        if col_idx == 1:
            matplotlib_support.turn_off_y_axis(axes)
        else:
            axes.set_ylabel('Mean normalized index')
        axes.set_ylim((0, 1))

        axes.legend(loc='lower right')
        axes.text(0.1, 0.9, structure, fontsize=12)
        matplotlib_support.set_axes_background(axes)
    fig.savefig(str(plot_path))
