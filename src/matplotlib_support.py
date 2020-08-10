
import numpy
import matplotlib
import math

import fig_style


def set_legend_background(legend, background_color='#ffffff', edge_color='#cccccc', alpha=0.7):
    frame = legend.get_frame()
    frame.set_color(background_color)
    frame.set_edgecolor(edge_color)
    frame.set_alpha(alpha)


def set_legend_marker_size(legend, marker_size):
    for handle in legend.legendHandles:
        handle._sizes = numpy.array([marker_size] * handle._sizes.size)


def turn_off_both_axis(axes):
    axes.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False,
                     left=False, right=False, labelleft=False,
                     )

def turn_off_x_axis(axes):
    axes.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)


def turn_off_y_axis(axes):
    axes.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)


def turn_off_both_axis(axes):
    turn_off_x_axis(axes)
    turn_off_y_axis(axes)
    axes.spines['left'].set_alpha(0)
    axes.spines['bottom'].set_alpha(0)


def set_y_ticks_right(axes):
    axes.tick_params(axis='y', which='both', right=True,
                     left=False, labelleft=False, labelright=True)
    axes.yaxis.set_label_position('right')


def set_axis_color(axes, color='black'):
    for child in axes.get_children():
        if isinstance(child, matplotlib.spines.Spine):
            child.set_color(color)


def set_axes_background(axes, color='white',
                        spine_left_color='grey', spine_bottom_color='grey'):
    axes.set_facecolor(color)
    axes.spines['left'].set_color(spine_left_color)
    axes.spines['bottom'].set_color(spine_bottom_color)


def turn_off_grid(axes):
    axes.grid(False)


def plot_legend(labels, colors, axes, fontize=fig_style.LEGEND_FONT_SIZE,
                nrows=None, marker='o',
                location='best'):

    legend_elements = []
    for label, color in zip(labels, colors):
        element = matplotlib.lines.Line2D([0], [0], marker=marker, color=color,
                         label=label,
                         markersize=fig_style.LEGEND_MARKER_SIZE, linestyle='')
        legend_elements.append(element)

    if nrows:
        ncol = len(legend_elements) // nrows
    else:
        ncol = 1
    legend = axes.legend(handles=legend_elements, prop={'size': fontize},
                         ncol=ncol, loc=location)
    set_legend_background(legend)


def set_y_ticks(tick_poss, tick_labels, axes, rotation=0, va='center', fontsize=10):
    axes.set_yticklabels(tick_labels, rotation=rotation, va=va, fontsize=fontsize)
    axes.set_yticks(tick_poss)


def set_x_ticks(tick_poss, tick_labels, axes, rotation=0, ha='right', fontsize=10):
    axes.set_xticklabels(tick_labels, rotation=rotation, ha=ha, fontsize=fontsize)
    axes.set_xticks(tick_poss)


def write_text_in_figure(text, x_pos, y_pos, fig, fontsize=12, zorder=200):
    axes = add_axes(fig, left_margin=0, right_margin=0,
                    top_margin=0, bottom_margin=0, zorder=zorder)
    axes.set_xlim((0, 1))
    axes.set_ylim((0, 1))
    axes.text(x_pos, y_pos, text, fontsize=fontsize)
    axes.set_facecolor('#00000000')
    axes.grid(False)


def add_axes(figure, row_idx=0, col_idx=0,
             left_margin=0.1, right_margin=0.02,
             top_margin=0.02, bottom_margin=0.07,
             axes_col_widths=None, axes_row_heights=None,
             projection=None, zorder=1):
    if axes_col_widths is None:
        axes_col_widths = [1]

    if axes_row_heights is None:
        axes_row_heights = [1]

    if isinstance(row_idx, int):
        row_idx = slice(row_idx, row_idx + 1)
    if isinstance(col_idx, int):
        col_idx = slice(col_idx, col_idx + 1)

    axes_col_widths = [width / sum(axes_col_widths) for width in axes_col_widths]
    axes_row_heights = [height / sum(axes_row_heights) for height in axes_row_heights]
    assert math.isclose(sum(axes_col_widths), 1)
    assert math.isclose(sum(axes_row_heights), 1)

    this_axes_width = sum(axes_col_widths[col_idx])
    prev_cols_width = sum(axes_col_widths[:col_idx.start])
    next_cols_width = sum(axes_col_widths[col_idx.stop:])

    this_axes_height = sum(axes_row_heights[row_idx])
    prev_rows_height = sum(axes_row_heights[:row_idx.start])
    next_rows_height = sum(axes_row_heights[row_idx.stop:])

    left_margin = left_margin * this_axes_width
    right_margin = right_margin * this_axes_width
    top_margin = top_margin * this_axes_height
    bottom_margin = bottom_margin * this_axes_height

    left = prev_cols_width + left_margin
    bottom = next_rows_height + bottom_margin
    width = this_axes_width - left_margin - right_margin
    height = this_axes_height - top_margin - bottom_margin
    return figure.add_axes((left, bottom, width, height), projection=projection, zorder=zorder)
