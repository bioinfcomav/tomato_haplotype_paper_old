
from itertools import combinations, count
from functools import partial
import unittest

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec

import colors
import matplotlib_support

INTERSECTION_AXES = 'intersection_axes'
SET_AXES = 'set_axes'
UPSET_AXES = 'upset_axes' 

ALL_INTERSECTIONS = 'all_intersections'
GIVEN_INTERSECTIONS = 'given_interesections'

LIGHT_GRAY = '#cccccc'
DARK_GRAY = '#333333'
DEFAULT_BAR_COLOR = 'tab:blue'
DEFAULT_EMPTY_DOT_SIZE = 50
DEFAULT_FULL_DOT_SIZE = 150
DEFAULT_AXIS_LABELS_FONT_SIZE = 15
DEFAULT_BARS_LEGEND_FONT_SIZE = int(DEFAULT_AXIS_LABELS_FONT_SIZE * 0.8)
UPSET_AXES_Y_AXIS_INCREMENT = 1.2
DEFAULT_COUNT_TYPE_COLOR_WHEEL = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                                  'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

def create_upset_axess(fig):
    gridspec = GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[0.8, 0.2])
    axess = {}
    axess[INTERSECTION_AXES] = fig.add_subplot(gridspec[0, 0])
    if False:
        axess[SET_AXES] = fig.add_subplot(gridspec[1, 0])
    axess[UPSET_AXES] = fig.add_subplot(gridspec[1, 0],
                                        sharex=axess[INTERSECTION_AXES])
                                        #sharey=axess[SET_AXES])
    return axess


def _sorting_key_function(intersection_data_labels, data_labels_order):
    num_data_labels = len(data_labels_order)
    numeric_index = [data_labels_order.index(data_label) + 1 for data_label in intersection_data_labels]
    numeric_index = [0] * (num_data_labels - len(numeric_index)) + numeric_index
    return numeric_index


def _get_uniq_data_labels(intersections_data_labels):
    return {data_label for intersection_data_labels in intersections_data_labels for data_label in intersection_data_labels}


def _get_sorted_intersections_data_labels(intersections_data_labels, intersections_to_plot, data_labels_order):

    uniq_data_labels = _get_uniq_data_labels(intersections_data_labels)

    given_intersections_data_labels = list(intersections_data_labels)
    if intersections_to_plot == GIVEN_INTERSECTIONS:
        intersections_data_labels = given_intersections_data_labels
    elif intersections_to_plot == ALL_INTERSECTIONS:
        items = uniq_data_labels
        intersections_data_labels = list([combination for num_items in range(1, len(items) + 1) for combination in combinations(items, num_items)])
    else:
        msg = 'ilegal value in intersections_to_plot'
        raise ValueError(msg)

    sorting_key_function = partial(_sorting_key_function, data_labels_order=data_labels_order)
    sorted_intersections_data_labels = sorted(intersections_data_labels, key=sorting_key_function)

    return sorted_intersections_data_labels


def _get_poss(num_items):
    return numpy.arange(num_items)


def plot_upset_dots(axes, present_intersections_data_labels,
                    sorted_intersections_data_labels,
                    sorted_uniq_data_labels,
                    empty_dot_size=DEFAULT_EMPTY_DOT_SIZE,
                    full_dot_size=DEFAULT_FULL_DOT_SIZE,
                    color_empty=LIGHT_GRAY,
                    color_full=DARK_GRAY):

    num_intersections = len(sorted_intersections_data_labels)
    num_data_labels = len(sorted_uniq_data_labels)

    x_poss = _get_poss(num_intersections)
    y_poss = _get_poss(num_data_labels)

    for intersection_idx, intersection_data_labels in enumerate(present_intersections_data_labels):
        x_pos = sorted_intersections_data_labels.index(intersection_data_labels)
        this_y_poss = [sorted_uniq_data_labels.index(data_label) for data_label in intersection_data_labels]
        if True:
            axes.scatter([x_pos] * len(this_y_poss), this_y_poss, s=full_dot_size,
                     color=color_full, edgecolors='none',
                     zorder=20)

    for y_pos in y_poss:
        if True:
            axes.scatter(x_poss, [y_pos] * num_intersections, s=empty_dot_size,
                     color=color_empty, edgecolors='none',
                     zorder=10)

        matplotlib_support.set_y_ticks(y_poss, sorted_uniq_data_labels, axes, fontsize=DEFAULT_AXIS_LABELS_FONT_SIZE)


def setup_axes(axess):
    if UPSET_AXES in axess:
        axes = axess[UPSET_AXES]
        axes.spines['left'].set_alpha(0)
        axes.spines['bottom'].set_alpha(0)
        axes.spines['right'].set_alpha(0)
        axes.spines['top'].set_alpha(0)
        axes.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        axes.tick_params(axis='y', which='both', right=False, left=False)
        bottom_lim, top_lim = axes.get_ylim()
        height = top_lim - bottom_lim
        middle = (top_lim + bottom_lim) / 2
        new_height = height * UPSET_AXES_Y_AXIS_INCREMENT
        bottom_lim = middle - new_height / 2
        top_lim = middle + new_height / 2
        axes.set_ylim(bottom_lim, top_lim)

    if INTERSECTION_AXES in axess:
        axes = axess[INTERSECTION_AXES]
        axes.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        axes.spines['right'].set_alpha(0)
        axes.spines['top'].set_alpha(0)


def plot_bar(axes, x_pos, height, color):
    axes.bar([x_pos], [height], bottom=None, color=color)


def plot_intersection_data(axes,
                           intersection_dict,
                           sorted_intersections_data_labels,
                           plot_intersection):

    for intersection_data_labels, values in intersection_dict.items():
        x_pos = sorted_intersections_data_labels.index(intersection_data_labels)
        plot_intersection(axes, x_pos, values)


def _plot_counts(intersection_dict, axess, plot_intersection, data_labels_order=None,
                 intersections_to_plot=GIVEN_INTERSECTIONS):

    if False:
        intersections_to_plot=ALL_INTERSECTIONS

    uniq_data_labels = _get_uniq_data_labels(intersection_dict.keys())

    if data_labels_order is None:
        data_labels_order = sorted(uniq_data_labels, key=str)

    sorted_uniq_data_labels = sorted(uniq_data_labels,
                                     key=lambda x: data_labels_order.index(x))

    intersections_data_labels = list(intersection_dict.keys())

    sorted_intersections_data_labels = _get_sorted_intersections_data_labels(intersections_data_labels,
                                                                             intersections_to_plot,
                                                                             data_labels_order)

    if UPSET_AXES in axess:
        upset_axes = axess[UPSET_AXES]
        plot_upset_dots(upset_axes,
                        intersections_data_labels,
                        sorted_intersections_data_labels,
                        sorted_uniq_data_labels)

    if INTERSECTION_AXES in axess:
        plot_intersection_data(axess[INTERSECTION_AXES],
                               intersection_dict,
                               sorted_intersections_data_labels,
                               plot_intersection)
    setup_axes(axess)


def plot_counts(intersection_dict, axess, data_labels_order=None,
                intersections_to_plot=GIVEN_INTERSECTIONS):
    plot_intersection = partial(plot_bar, color=DEFAULT_BAR_COLOR)
    _plot_counts(intersection_dict=intersection_dict,
                 axess=axess, data_labels_order=data_labels_order,
                 intersections_to_plot=intersections_to_plot,
                 plot_intersection=plot_intersection)


def plot_stacked_bar(axes, x_pos, heights_by_type,
                     count_type_color_schema, sorted_count_types):

    bottom = None
    for count_type in sorted_count_types:
        if count_type not in heights_by_type:
            continue
        height = heights_by_type[count_type]
        color = count_type_color_schema[count_type]
        axes.bar([x_pos], [height], bottom=bottom, color=color)
        if bottom is None:
            bottom = height
        else:
            bottom += height


def plot_counts_per_type(intersection_dict, axess, data_labels_order=None,
                         intersections_to_plot=GIVEN_INTERSECTIONS,
                         sorted_count_types=None,
                         count_type_color_schema=None):

    if sorted_count_types is None:
        count_types = {count_type for counts_per_type in intersection_dict.values() for count_type in counts_per_type.keys()}
        sorted_count_types = sorted(count_types, key=str)

    if count_type_color_schema is None:
        count_type_color_schema = colors.ColorSchema(wheel=DEFAULT_COUNT_TYPE_COLOR_WHEEL)

    plot_intersection = partial(plot_stacked_bar, count_type_color_schema=count_type_color_schema,
                                sorted_count_types=sorted_count_types)

    _plot_counts(intersection_dict=intersection_dict,
                 axess=axess, data_labels_order=data_labels_order,
                 intersections_to_plot=intersections_to_plot,
                 plot_intersection=plot_intersection)
    
    matplotlib_support.plot_legend(sorted_count_types,
                                   [count_type_color_schema[count_type] for count_type in sorted_count_types],
                                   axess[INTERSECTION_AXES],
                                   fontize=DEFAULT_BARS_LEGEND_FONT_SIZE)
    return {'axess': axess}


class PlotCountsTest(unittest.TestCase):
    def test_plot_counts(self):

        fig = Figure()
        FigureCanvas(fig)

        axess = create_upset_axess(fig)

        intersection_dict = {('a', ): 1,
                             ('c',): 5,
                             ('a', 'b'): 3,
                             ('c', 'b'): 4}
        plot_counts(intersection_dict, axess)
        fig.savefig('../tmp/upset.svg')

    def test_plot_counts_per_type(self):

        fig = Figure()
        FigureCanvas(fig)

        axess = create_upset_axess(fig)

        intersection_dict = {('a', ): {'type1': 1},
                             ('c',): {'type2': 5},
                             ('a', 'b'): {'type1': 1, 'type2': 2},
                             ('c', 'b'): {'type1': 2, 'type2': 2}}
        plot_counts_per_type(intersection_dict, axess)
        fig.savefig('../tmp/upset2.svg')


if __name__ == '__main__':
    unittest.main()
