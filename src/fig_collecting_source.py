
import config

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import morphological
import colors
import matplotlib_support


if __name__ == '__main__':
    original_data = morphological.read_morphological_data()

    plot_path = config.FIG_COLLECTING_SOURCES

    source_order = ['natural', 'natural_or_disturbed', 'disturbed',
                    'weed_or_semi-cultivated', 'semi-cultivated', 'semi-cultivated_or_cultivated',
                    'cultivated', 'market']

    taxon_mapping = {('SLC', 'CRI'): ('SLC', 'MA'),
                     ('SLC', 'MEX'): ('SLC', 'MA')}
    taxon_order = (('SP', 'PER'), ('SP', 'ECU'), ('SLC', 'COL'),
                   ('SLC', 'CRI'), ('SLC', 'MEX'), ('SLC', 'MA'), ('SLC', 'PER'),
                   ('SLC', 'PER_N'), ('SLC', 'ECU'), ('SLL', 'MEX'))

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes_row_heights = [0.2, 0.8]
    axes = matplotlib_support.add_axes(fig, row_idx=1, col_idx=0, axes_row_heights=axes_row_heights,
                                       bottom_margin=0.2)
    res = morphological.plot_collecting_sources(original_data, axes, color_schema=colors.ColorSchema(colors.SOURCE_COLORS),
                                                source_order=source_order,
                                                taxon_mapping=taxon_mapping,
                                                taxon_order=taxon_order)
    matplotlib_support.set_axes_background(axes)

    axes = matplotlib_support.add_axes(fig, row_idx=0, col_idx=0, axes_row_heights=axes_row_heights)
    axes.legend(res['artists'], res['labels'], ncol=3)
    matplotlib_support.set_axes_background(axes)
    matplotlib_support.turn_off_both_axis(axes)
    fig.savefig(plot_path)
