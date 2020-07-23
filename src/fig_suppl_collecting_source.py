
import config

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import morphological
import colors
import matplotlib_support

if __name__ == '__main__':
    original_data = morphological.read_morphological_data()

    out_dir = config.FIGURES_DIR
    plot_path = out_dir / 'collecting_sources.svg'

    source_order = ['natural', 'natural_or_disturbed', 'disturbed',
                    'weed_or_semi-cultivated', 'semi-cultivated', 'semi-cultivated_or_cultivated',
                    'cultivated', 'market']

    taxon_mapping = {('SLC', 'CRI'): ('SLC', 'MA'),
                     ('SLC', 'MEX'): ('SLC', 'MA')}
    taxon_order = ('SP', 'PER'), ('SP', 'ECU'), ('SLC', 'COL'), ('SLC', 'CRI'), ('SLC', 'MEX'), ('SLC', 'MA'), ('SLC', 'PER'), ('SLC', 'ECU'), ('SLL', 'MEX')

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)
    morphological.plot_collecting_sources(original_data, axes, color_schema=colors.ColorSchema(colors.SOURCE_COLORS),
                            source_order=source_order, taxon_mapping=taxon_mapping, taxon_order=taxon_order)
    matplotlib_support.set_axes_background(axes)
    axes.legend()
    fig.tight_layout()
    fig.savefig(plot_path)
