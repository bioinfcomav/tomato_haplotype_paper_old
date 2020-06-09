
import config

from collections import defaultdict
from pprint import pprint

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.image import imread

import cartopy.crs as ccrs

import passport
from colors import (CLASSIFICATION_RANK1_COLORS, CLASSIFICATION_RANK2_COLORS,
                    ColorSchema, modify_rgb_hex_color_hsv,
                    modify_rgb_image_hsv, rgb_image_to_rgb_0_1)
from labels import get_long_label, LABELS
from plot import get_sorted_legend_handles


def plot_geo_map(samples, axes, classifications=None, color_schema=None,
                 plot_sample_ids=False,
                 longitude_range=None, latitude_range=None,
                 images=None, draw_coastlines=False):

    if images is None:
        images = []


    if classifications is None:
        classifications = {}

    if color_schema is None:
        color_schema = ColorSchema()

    if draw_coastlines:
        axes.coastlines(zorder=1)

    axes.background_patch.set_visible(False)

    samples_to_plot = defaultdict(dict)
    for sample_id, sample_info in samples.items():
        if 'latitude' not in sample_info:
            continue

        classification = classifications.get(sample_id)
    
        if classification not in samples_to_plot:
            samples_to_plot[classification] = {'samples': [], 'latitudes': [], 'longitudes': []}
        longitude = sample_info['longitude']
        latitude = sample_info['latitude']

        if longitude is None or latitude is None:
            continue

        if longitude_range is not None:
            if longitude < longitude_range[0] or longitude > longitude_range[1]:
                continue

        if latitude_range is not None:
            if latitude < latitude_range[0] or latitude > latitude_range[1]:
                continue
    
        samples_to_plot[classification]['samples'].append(sample_id)
        samples_to_plot[classification]['latitudes'].append(latitude)
        samples_to_plot[classification]['longitudes'].append(longitude)

    order = {pop: idx for idx, pop in enumerate(reversed(LABELS.keys()))}
    sorted_pops = sorted(samples_to_plot.keys(), key=lambda x: order.get(x, 1000))

    for pop in sorted_pops:
        samples_info = samples_to_plot[pop]
        color = color_schema[pop]
        #color = modify_rgb_hex_color_hsv(color, luminosity_addition=0.2)
        border_color = modify_rgb_hex_color_hsv(color, luminosity_addition=-0.2)

        label = get_long_label(pop)

        axes.scatter(samples_info['longitudes'], samples_info['latitudes'],
                     label=label, c=color, zorder=10, s=120,
                     edgecolors=border_color)

        if plot_sample_ids:
            for x, y, text in zip(samples_info['longitudes'],
                                  samples_info['latitudes'],
                                  samples_info['samples']):
                axes.text(x, y, text)

    x_lims = axes.get_xlim()
    y_lims = axes.get_ylim()

    image_proj = ccrs.PlateCarree()
    for image in images:
        if image.get('ignore'):
            continue
        rgb_image = imread(str(image['path']))
        rgb_image = rgb_image_to_rgb_0_1(rgb_image)
        if 'hsv_modification' in image:
            rgb_image = modify_rgb_image_hsv(rgb_image, **image['hsv_modification'])
        axes.imshow(rgb_image,
                    origin='upper', transform=image_proj, 
                    extent=image['extent'],
                    zorder=image['zorder'],
                    )

    axes.set_xlim(x_lims)
    axes.set_ylim(y_lims)

    legend_handles_and_labels = get_sorted_legend_handles(axes)

    return {'legend_handles_and_labels': legend_handles_and_labels}


def plot_geo_rank1_for_main_pops(passports, axes=None, plot_legend=True):

    colors = ColorSchema(CLASSIFICATION_RANK1_COLORS)

    rank = 'rank1'

    revelant_pops = {'rank1': ['sp_ec', 'sll_mx', 'slc_ec', 'sp_pe_inter-andean',
                    'slc_ma', 'slc_pe', 'sp_pe', 'slc_co']}
    revelant_pops = revelant_pops[rank]

    passports = {sample_id: samples_info for sample_id, samples_info in passports.items() if samples_info.get('classification', {}).get(rank) in revelant_pops}

    classifications = {sample_id: passport.get('classification', {}).get(rank) for sample_id, passport in passports.items()}

    hypothesis_path = config.HYPOTHESIS_PNG

    if axes is None:
        plot_path = config.GEOGRAPHIC_FIGURE_DIR / f'geo_map.svg'
        fig = Figure((10, 15))
        FigureCanvas(fig) # Don't remove it or savefig will fail later
        axes = fig.add_subplot(111, projection=ccrs.PlateCarree(), zorder=1)
        savefig = True
    else:
        savefig = False

    res = plot_geo_map(passports, axes=axes,
                       classifications=classifications,
                       color_schema=colors,
                       plot_sample_ids=False, longitude_range=(-116, -60),
                       images=[{'ignore': False,
                               'path': config.NE_BACKGROUND_CUT_PNG,
                               'extent': (-111.5, -66.9, -19.2, 28.2),
                               'zorder': 1,
                               'hsv_modification': {'luminosity_addition': 0.2,
                                                    'saturation_addition': -0.2}
                               },
                               {'path': hypothesis_path,
                               'extent': (-106.31449895265104, -69.52907660642875, -19.596286239604375, 27.67473222177241),
                               'zorder': 20}])

    if plot_legend:
        legend_handles, legend_labels = res['legend_handles_and_labels']
        handle_and_labels = zip(legend_handles, legend_labels)

        order = {label: idx for idx, label in enumerate(LABELS.values())}

        handle_and_labels = sorted(handle_and_labels, key=lambda t: t[1])
        handle_and_labels = sorted(handle_and_labels, key=lambda t: order.get(t[1], 1000))
        handles, labels = zip(*handle_and_labels)
        axes.legend(handles, labels,
                    prop={'size': 17},
                    loc='lower left')

    if savefig:
        fig.savefig(str(plot_path))


def plot_geo_supplemental_rank2_for_all_pops(passports):

    colors = ColorSchema(CLASSIFICATION_RANK2_COLORS)

    rank = 'rank2'

    pops_to_exclude = ['slc_co', 'sp_x_sl_cherry_cultivars', 'sll_vint', 'sll_vint_small',
                       'sll_moneymaker_ailsacraig', 'sll_oldbrooks',  'sll_oldgerman',
                       'sll_modern']
    passports = {sample_id: samples_info for sample_id, samples_info in passports.items() if samples_info.get('classification', {}).get('rank2') not in pops_to_exclude}

    classifications = {sample_id: passport.get('classification', {}).get(rank) for sample_id, passport in passports.items()}

    fig = Figure((10, 20))
    FigureCanvas(fig) # Don't remove it or savefig will fail later

    grid_spec = fig.add_gridspec(nrows=2, ncols=2,
                                 height_ratios=(3.1, 1),
                                 hspace=0.05, wspace=0.03)

    mesoamerican_axes = fig.add_subplot(grid_spec[0, 0], projection=ccrs.PlateCarree(), zorder=1)

    plot_path = config.GEOGRAPHIC_FIGURE_DIR / f'geo_map_supplemental.svg'

    plot_background = False
    draw_coastlines = False

    res = plot_geo_map(passports, axes=mesoamerican_axes,
                       classifications=classifications,
                       color_schema=colors,
                       draw_coastlines=draw_coastlines,
                       longitude_range=(-112, -81),
                       latitude_range=(9, 26),
                       images=[{'ignore': False,
                                'path': config.NE_BACKGROUND_CUT_PNG,
                                'extent': (-111.5, -66.9, -19.2, 28.2),
                                'zorder': 1,
                                'hsv_modification': {'luminosity_addition': 0.2,
                                                     'saturation_addition': -0.2}}]
                       )
    mesoamerican_axes.text(-110, 9, 'A', {'size': 25})

    legend_handles1, legend_labels1 = res['legend_handles_and_labels']

    andean_axes = fig.add_subplot(grid_spec[:, 1], projection=ccrs.PlateCarree(), zorder=1)
    andean_axes.text(-72, 0.5, 'B', {'size': 25})

    res = plot_geo_map(passports, axes=andean_axes,
                       classifications=classifications,
                       color_schema=colors,
                       draw_coastlines=draw_coastlines,
                       longitude_range=(-81, -70),
                       latitude_range=(-17, 1.5),
                       images=[{'ignore': False,
                                'path': config.NE_BACKGROUND_CUT_PNG,
                                'extent': (-111.5, -66.9, -19.2, 28.2),
                                'zorder': 1,
                                'hsv_modification': {'luminosity_addition': 0.2,
                                                     'saturation_addition': -0.2}}]
                       )
    legend_handles2, legend_labels2 = res['legend_handles_and_labels']

    handles_and_labels = zip(legend_handles1 + legend_handles2,
                                legend_labels1 + legend_labels2)
    labels_seen = set()
    handle_and_labels = []
    for handle_and_label in handles_and_labels:
        if handle_and_label[1] in labels_seen:
            continue
        labels_seen.add(handle_and_label[1])
        handle_and_labels.append(handle_and_label)

    order = {label: idx for idx, label in enumerate(LABELS.values())}

    handle_and_labels = sorted(handle_and_labels, key=lambda t: t[1])
    handle_and_labels = sorted(handle_and_labels, key=lambda t: order.get(t[1], 1000))
    handles, labels = zip(*handle_and_labels)

    #legend_axes = fig.add_subplot(grid_spec[1, 0])

    andean_axes.legend(handles, labels, prop={'size': 9.4},
                       bbox_to_anchor=(-0.01, 0.61, 0, 0),
                       ncol=2)

    fig.savefig(str(plot_path))


if __name__ == '__main__':
    samples = passport.get_sample_passports()

    if True:
        plot_geo_rank1_for_main_pops(samples)
    if True:
        plot_geo_supplemental_rank2_for_all_pops(samples)
