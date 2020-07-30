
import config

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import cartopy.crs as ccrs

import morphological
import passport
import labels
import colors
import matplotlib_support
import plot_geo_map
import morphological_progression


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
                                       right_margin=0, left_margin=0, bottom_margin=0.05,
                                       projection=ccrs.PlateCarree())

    passports = {acc: {'longitude': data.get('location', {}).get('longitude'), 'latitude': data.get('location', {}).get('latitude')} for acc, data in original_data.items()}
    plot_geo_map.plot_geo_map(passports, axes=axes,
                              classifications=morpho_classification,
                              color_schema=color_schema,
                              plot_sample_ids=False, longitude_range=(-116, -60),
                              images=[{'ignore': False,
                                       'path': config.NE_BACKGROUND_CUT_PNG,
                                       'extent': (-111.5, -65.0, -19.5, 30),
                                       'zorder': 1,
                                       'hsv_modification': {'luminosity_addition': 0.2,
                                                            'saturation_addition': -0.2}}])
    axes.legend()

    morpho_classes = ['sp_pe', 'sp_pe_inter-andean', 'sp_ec',
                      'slc_ec', 'slc_ma', 'slc_pe', 'sll']
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
