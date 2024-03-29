
import config

from collections import defaultdict

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import morphological
import matplotlib_support
import colors


def plot_morphological_progression(morpho_data, classes_for_accessions, axes,
                                   sorted_classes=None, traits=None,
                                   normalize=True,
                                   background_colors=None,
                                   bar_alpha=1, linestyle='-', marker='',
                                   color_schema=None,
                                   add_label=True):

    if color_schema is None:
        color_schema = colors.ColorSchema()

    if sorted_classes is None:
        sorted_classes = sorted({klass for klass in classes_for_accessions.values()}, key=str)

    accs_by_morpho_class = defaultdict(list)
    for acc, morpho_class in classes_for_accessions.items():
        accs_by_morpho_class[morpho_class].append(acc)

    means_per_morpho_class = {}
    for morpho_class, accs in accs_by_morpho_class.items():
        common_accs = set(accs).intersection(morpho_data.index)
        morpho_data_for_this_morpho_class = morpho_data.reindex(common_accs)
        means = morpho_data_for_this_morpho_class.mean()
        for trait, mean in means.to_dict().items():
            if trait not in means_per_morpho_class:
                means_per_morpho_class[trait] = {}
            means_per_morpho_class[trait][morpho_class] = mean

    if traits is None:
        traits = list(means_per_morpho_class.keys())

    edge_poss = numpy.arange(len(sorted_classes) + 1)
    x_values = (edge_poss[:-1] + edge_poss[1:]) / 2

    for trait in traits:
        y_values = numpy.array([means_per_morpho_class[trait][morpho_class] for morpho_class in sorted_classes])
        y_values[numpy.isnan(y_values)] = 0
        if normalize:
            y_values = numpy.array(y_values)
            min_ = numpy.min(y_values)
            max_ = numpy.max(y_values)
            y_values = (y_values -min_) / (max_ - min_)

        color = color_schema[trait]

        if add_label:
            axes.plot(x_values, y_values, label=trait, zorder=10,
                    linestyle=linestyle, marker=marker, color=color)
        else:
            axes.plot(x_values, y_values, zorder=10,
                    linestyle=linestyle, marker=marker, color=color)

    if background_colors:
        y_lims = axes.get_ylim()
        height = y_lims[1], y_lims[0]
        width = edge_poss[1] - edge_poss[0]
        bottom = y_lims[0]
        for morpho_class, x0 in zip(sorted_classes, edge_poss[:-1]):
            color = background_colors[morpho_class]
            axes.bar([x0], height, width=width, bottom=bottom, align='edge', zorder=1, color=color, alpha=bar_alpha)
        axes.set_xlim((edge_poss[0], edge_poss[-1]))

    matplotlib_support.set_x_ticks(x_values, sorted_classes, axes, rotation=45)
    return {'color_schema': color_schema}


if __name__ == '__main__':
    #original_data = morphological.read_morphological_data()
    data = morphological.get_morphological_table_for_ordinal_traits()

    morpho_classification = morphological.read_morphological_classification()

    morpho_classes = ['sp_peru', 'sp_inter_andean', 'sp_ecu',
                      'slc_ecu', 'slc_ma', 'slc_pe', 'sll']

    out_dir = config.MORPHOLOGICAL_DIR

    plot_path = out_dir / 'morphological_progression_all.svg'
    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)
    plot_morphological_progression(data, morpho_classification, axes, sorted_morpho_classes=morpho_classes)
    axes.legend()
    fig.tight_layout()
    fig.savefig(str(plot_path))

    plot_path = out_dir / 'morphological_progression.svg'
    fig = Figure((5, 10))
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes1 = fig.add_subplot(411)
    traits = ['presence_of_irregular_inflorescence', 'ribbing', 'style_exsertion',
              'petal_position', 'stem_width', 'stem_hairiness', 'leaf_type',
              'petal_width', 'fruit_size', 'leaflet_margin', 'inflorescence_ordinal']

    traits = ['presence_of_irregular_inflorescence', 'style_exsertion', 'style_curvature',
              'petal_position', 'petal_width', 'inflorescence_ordinal']
    plot_morphological_progression(data, morpho_classification, axes1, sorted_morpho_classes=morpho_classes,
                                   traits=traits)
    axes2 = fig.add_subplot(412)
    traits = ['stem_width', 'stem_hairiness']
    plot_morphological_progression(data, morpho_classification, axes2, sorted_morpho_classes=morpho_classes,
                                   traits=traits)

    axes3 = fig.add_subplot(413)
    traits = ['leaf_type', 'leaflet_margin']
    plot_morphological_progression(data, morpho_classification, axes3, sorted_morpho_classes=morpho_classes,
                                   traits=traits)

    axes4 = fig.add_subplot(414)
    traits = ['ribbing', 'fruit_size', 'peanut_fruit_shape', 'stripy_fruit']
    plot_morphological_progression(data, morpho_classification, axes4, sorted_morpho_classes=morpho_classes,
                                   traits=traits)

    axes1.legend()
    axes2.legend()
    axes3.legend()
    axes4.legend()
    matplotlib_support.turn_off_x_axis(axes1)
    matplotlib_support.turn_off_x_axis(axes2)
    matplotlib_support.turn_off_x_axis(axes3)
    matplotlib_support.set_axes_background(axes1)
    matplotlib_support.set_axes_background(axes2)
    matplotlib_support.set_axes_background(axes3)
    matplotlib_support.set_axes_background(axes4)
    axes1.set_ylabel('Mean normalized index')
    axes2.set_ylabel('Mean normalized index')
    axes3.set_ylabel('Mean normalized index')
    axes4.set_ylabel('Mean normalized index')
    fig.tight_layout()
    fig.savefig(str(plot_path))
