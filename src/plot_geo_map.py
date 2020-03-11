
import config

from collections import defaultdict
from pprint import pprint

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.image import imread

import cartopy.crs as ccrs

import passport
from colors import (CLASSIFICATION_COLORS, ColorSchema, modify_rgb_hex_color_hsv,
                    modify_rgb_image_hsv, rgb_image_to_rgb_0_1)
from labels import get_long_label
from plot import get_sorted_legend_handles


def plot_geo_map(samples, classification_rank, axes, plot_sample_ids=False,
                 longitude_range=None, latitude_range=None,
                 images=None, draw_coastlines=False, plot_legend=False):

    if images is None:
        images = []

    if draw_coastlines:
        axes.coastlines(zorder=1)

    axes.background_patch.set_visible(False)

    colors = ColorSchema(CLASSIFICATION_COLORS)

    samples_to_plot = defaultdict(dict)
    for sample_id, sample_info in samples.items():
        if 'latitude' not in sample_info:
            continue

        classification = str(sample_info.get('classification', {}).get(classification_rank))
        if classification not in samples_to_plot:
            samples_to_plot[classification] = {'samples': [], 'latitudes': [], 'longitudes': []}

        longitude = sample_info['longitude']
        latitude = sample_info['latitude']

        if longitude_range is not None:
            if longitude < longitude_range[0] or longitude > longitude_range[1]:
                continue

        if latitude_range is not None:
            if latitude < latitude_range[0] or latitude > latitude_range[1]:
                continue
    
        samples_to_plot[classification]['samples'].append(sample_id)
        samples_to_plot[classification]['latitudes'].append(latitude)
        samples_to_plot[classification]['longitudes'].append(longitude)

    for classification, samples_info in samples_to_plot.items():
        color = colors[classification]
        #color = modify_rgb_hex_color_hsv(color, luminosity_addition=0.2)
        border_color = modify_rgb_hex_color_hsv(color, luminosity_addition=-0.2)

        label = get_long_label(classification)

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
        print(image)
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

    if plot_legend:
        axes.legend(*get_sorted_legend_handles(axes),
                    prop={'size': 17},
                    loc='lower left')


def plot_geo_rank1_for_main_pops(samples):

    rank = 'rank1'

    revelant_pops = {'rank1': ['sp_ec', 'sll_mx', 'slc_ec', 'sp_pe_inter-andean',
                    'slc_ma', 'slc_pe', 'sp_pe', 'slc_co']}
    revelant_pops = revelant_pops[rank]

    passports = {sample_id: samples_info for sample_id, samples_info in samples.items() if samples_info.get('classification', {}).get(rank) in revelant_pops}

    plot_path = config.GEOGRAPHIC_FIGURE_DIR / f'geo_map.svg'
    hypothesis_path = config.HYPOTHESIS_PNG

    fig = Figure((10, 15))
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111, projection=ccrs.PlateCarree(), zorder=1)

    plot_geo_map(passports, axes=axes,
                classification_rank=rank,
                plot_sample_ids=False, longitude_range=(-116, -60),
                images=[{#'ignore': True,
                        'path': config.NE_BACKGROUND_TIFF,
                        'extent': [-180, 180, -90, 90],
                        'zorder': 1,
                        'hsv_modification': {'luminosity_addition': 0.2,
                                            'saturation_addition': -0.2}
                        },
                        {'path': hypothesis_path,
                        'extent': (-106.31449895265104, -69.52907660642875, -19.596286239604375, 27.67473222177241),
                        'zorder': 20}])
    fig.tight_layout()
    fig.savefig(str(plot_path))


def plot_geo_supplemental_rank2_for_all_pops(passports):

    fig = Figure((7, 15))
    FigureCanvas(fig) # Don't remove it or savefig will fail later

    grid_spec = fig.add_gridspec(nrows=2, ncols=1,
                                 width_ratios=(1,), height_ratios=(1, 2.6),
                                 hspace=0.05)

    mesoamerican_axes = fig.add_subplot(grid_spec[0, 0], projection=ccrs.PlateCarree(), zorder=1)

    plot_path = config.GEOGRAPHIC_FIGURE_DIR / f'geo_map_supplemental.svg'
    rank = 'rank2'

    plot_background = True

    plot_geo_map(passports, classification_rank=rank, axes=mesoamerican_axes,
                draw_coastlines=True,
                longitude_range=(-112, -81),
                latitude_range=(9, 26),
                plot_legend=False,
                images=[{'ignore': not(plot_background),
                         'path': config.NE_BACKGROUND_TIFF,
                         'extent': [-180, 180, -90, 90],
                         'zorder': 1,
                         'hsv_modification': {'luminosity_addition': 0.2,
                                             'saturation_addition': -0.2}}],
                )

    andean_axes = fig.add_subplot(grid_spec[1, 0], projection=ccrs.PlateCarree(), zorder=1)

    plot_geo_map(passports, classification_rank=rank, axes=andean_axes,
                 draw_coastlines=True,
                 longitude_range=(-81, -70),
                 latitude_range=(-17, 1.5),
                 plot_legend=False,
                 images=[{'ignore': not(plot_background),
                         'path': config.NE_BACKGROUND_TIFF,
                         'extent': [-180, 180, -90, 90],
                         'zorder': 1,
                         'hsv_modification': {'luminosity_addition': 0.2,
                                             'saturation_addition': -0.2}}]
                 )


    fig.tight_layout()
    fig.savefig(str(plot_path))


if __name__ == '__main__':
    samples = passport.get_sample_passports()

    if False:
        plot_geo_rank1_for_main_pops(samples)
    else:
        plot_geo_supplemental_rank2_for_all_pops(samples)
