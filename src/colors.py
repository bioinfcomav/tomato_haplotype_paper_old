
import numpy

from matplotlib.colors import LinearSegmentedColormap, rgb_to_hsv, hsv_to_rgb
{
                         'slc_ma': '#20B2AA', # lightseagreen
                         'sp_ec': '#808000', # olive
                         'sp_pe': '#556B2F',   # darkolive
                         'slc_co': '#2bae2b',
                         'sp_pe_inter-andean': '#36835e', 

                         'sll_mx': '#FF6347', #tomato
                        'slc_ec': '#ffabb1',
                        'slc_pe': '#dca675'
}

CLASSIFICATION_RANK1_COLORS = {'sll_mex': '#ff0000',
                         'slc_ca_small': '#FF69B4',
                         'slc_tarapoto': '#57acff',
                         'slc_moyobamba': '#0043ff',
                         'slc_ecu_center': '#3CB371', 
                         'slc_ecu_north': '#026b17',
                         'slc_peru_south': '#fffa05',
                         'sp_peru': '#cccccc',
                         'sp_montane': '#888888',
                         'sp_ecu': '#a35c11',
                         'gal': '#FFA500',
                         'sll_vint': '#CD5C5C',
                         'sll_vint_small': '#E9967A',
                        
                         'slc_ma': '#b3b300', 
                         'sp_ec': '#2bae2b',
                         'sp_pe': '#556B2F',   # darkolive
                         'slc_co': '#f5ea6f',
                         'sp_pe_inter-andean': '#20B2AA', 

                         'sll_mx': '#FF6347', #tomato
                         'slc_ec': '#ffc087',
                         'slc_pe': '#ffabb1'

                         }

CLASSIFICATION_RANK2_COLORS = {
'sp_pe_desert': '#556B2F',
'sp_pe_n_hills': '#8FBC8F',
'sp_pe_n_inter-andean': '#20B2AA',

'sp_ec_s_dry_forest': '#9ACD32',
'sp_ec_n_wet_forest': '#2bae2b',

'sp_pe_desert_x_sp_pe_n_inter_andean': '#D8D8D8',
'sp_pe_x_sp_ec': '#D8D8D8',

'slc_ec_n_600m': '#D2B48C',
'slc_ec_c_800m': '#A46D5A',
'slc_ec_s_1000m': '#7A4C1E',
'slc_ec_guayaquil': '#b83434',

#'slc_co': '#f5ea6f',
'slc_mx_sinaloa': '#ffa500',
'slc_ca': '#ffee58',
'slc_mx': '#b3b300',

'sp_x_sl': '#D8D8D8',
'sp_x_sl_cherry_cultivars': '#D8D8D8',

'slc_pe_n_400m': '#ffabb1',
'slc_pe_n_1000m': '#FF69B4',
'slc_pe_c': '#DDA0DD',

'sll_mx': '#FF6347',
'slc_world': '#DC143C',

#'sll_vint': '',
#'sll_vint_small': '',
#'sll_moneymaker_ailsacraig': '',
#'sll_oldbrooks': '',
#'sll_oldgerman': '',
#'sll_modern': ''
}

SMALL_COLOR_WHEEL = ['#57acff', '#ff6347', '#3CB371']
COLOR_WHEEL = ['#ff6347', '#57acff', '#8FBC8F', '#FF69B4', '#CD5C5C',
               '#80ced6', '#808000', '#ba68c8', '#FFA07A', '#92a8d1',
               '#3CB371', '#D2B48C', '#667292']

color_to_1_0 = lambda x: x / 255
hex_to_rgb_1_0 = lambda h: tuple(int(h[i:i+2], 16) / 255 for i in (1, 3, 5))

white = 255, 255, 255
black = 0, 0, 0
palepink =  253, 237, 236
salmon2 = 250, 219, 216
salmon = 241, 148, 138
salmon3 = 245, 170, 160
paleorange =  240, 178, 122
greenblue = 23, 165, 137
darkgreen = 11, 83, 69
color2 = 236, 231, 242
color3 = 208, 209, 230
color4 = 166, 189, 219
color5 = 116, 169, 207
color6 = 54, 144, 192
color7 = 5, 112, 176
color8 = 4, 90, 141
color9 = 2, 56, 88

white = list(map(color_to_1_0, white))
black = list(map(color_to_1_0, black))
palepink = list(map(color_to_1_0, palepink))
color2 = list(map(color_to_1_0, color2))
color3 = list(map(color_to_1_0, color3))
color4 = list(map(color_to_1_0, color4))
color5 = list(map(color_to_1_0, color5))
color6 = list(map(color_to_1_0, color6))
color7 = list(map(color_to_1_0, color7))
color8 = list(map(color_to_1_0, color9))
color9 = list(map(color_to_1_0, color8))
gwas_color_list = [(0, palepink), (0.125, color2), (0.25, color3), (0.375, color4), (0.5, color5), (0.625, color6), (0.75, color7), (0.875, color8), (1, color9)]
gwas_color_list_r = [(0, color9), (0.125, color8), (0.25, color7), (0.375, color6), (0.5, color5), (0.625, color4), (0.75, color3), (0.875, color2), (1, palepink)]
GWAS_CMAP = LinearSegmentedColormap.from_list('rg', gwas_color_list, N=256)

color1 = hex_to_rgb_1_0('#1c0c74')
color2 = hex_to_rgb_1_0('#2979bb')
color3 = hex_to_rgb_1_0('#8cb4c1')
color4 = hex_to_rgb_1_0('#d8cbe5')
color5 = hex_to_rgb_1_0('#fbf4f4')
pink_blue_list = [(0, color1), (0.2, color2), (0.4, color3), (0.6, color4), (1.0, color5)]
pink_blue_list_r = [(0, color5), (0.25, color4), (0.5, color3), (0.75, color2), (1.0, color1)]
pink_blue_list_r = [(0, white),
                    (0.05, color5),
                    (0.15, color4),
                    (0.35, color3),
                    (0.55, color2),
                    (0.75, color1),
                    (1.0, black)]
PINK_BLUE_CMAP = LinearSegmentedColormap.from_list('rg',
                                                    pink_blue_list,
                                                    N=256)
PINK_BLUE_CMAP_R = LinearSegmentedColormap.from_list('rg',
                                                    pink_blue_list_r,
                                                    N=256)

class ColorSchema:
    def __init__(self, colors_used=None, wheel=COLOR_WHEEL):
        self._wheel = COLOR_WHEEL
        if colors_used is None:
            colors_used = {}
        self.colors_used = colors_used.copy()

    def __getitem__(self, key):
        if key in self.colors_used:
            return self.colors_used[key]

        wheel = self._wheel
        color_idx = (len(self.colors_used) + 1) % len(wheel)
        color = wheel[color_idx]
        self.colors_used[key] = color
        return color


def hex_rgb_color_to_rgb_tuple(color_hex_rgb):
    return tuple(int(color_hex_rgb.lstrip('#')[i:i+2], 16) / 255 for i in (0, 2, 4))


def hex_rgb_tuple_to_rgb_hex(color_tuple_rgb):
    return '#%02x%02x%02x' % (int(color_tuple_rgb[0] * 255), int(color_tuple_rgb[1] * 255), int(color_tuple_rgb[2] * 255))


def lower_color_luminosity(color_hex_rgb, lower_rate=0.1):
    rgb_color_tuple = hex_rgb_color_to_rgb_tuple(color_hex_rgb)
    hsv_color = rgb_to_hsv(rgb_color_tuple)

    lowered_luminosity = hsv_color[2] - lower_rate
    if lowered_luminosity < 0:
        lowered_luminosity = 0
    hsv_color = hsv_color[0], hsv_color[1], lowered_luminosity

    rgb_color = hsv_to_rgb(hsv_color)
    return hex_rgb_tuple_to_rgb_hex(rgb_color)


def modify_rgb_hex_color_hsv(color_hex_rgb, luminosity_addition=0., hue_addition=0., saturation_addition=0.):
    rgb_color_tuple = hex_rgb_color_to_rgb_tuple(color_hex_rgb)
    hsv_color = rgb_to_hsv(rgb_color_tuple)

    hue = hsv_color[0] + hue_addition
    saturation = hsv_color[1] + saturation_addition
    luminosity = hsv_color[2] + luminosity_addition

    hsv_color = numpy.array((hue, saturation, luminosity))

    hsv_color[hsv_color < 0] = 0
    hsv_color[hsv_color > 1] = 1

    rgb_color = hsv_to_rgb(hsv_color)
    return hex_rgb_tuple_to_rgb_hex(rgb_color)


def modify_rgb_image_hsv(image_array, luminosity_addition=0., hue_addition=0., saturation_addition=0.):
    img_shape = image_array.shape[:2]

    hue_alteration = numpy.full((img_shape[0], img_shape[1]), hue_addition)
    saturation_alteration = numpy.full((img_shape[0], img_shape[1]), saturation_addition)
    luminosity_alteration = numpy.full((img_shape[0], img_shape[1]), luminosity_addition)
    alteration = numpy.stack([hue_alteration, saturation_alteration, luminosity_alteration], axis=2)

    hsv_image = rgb_to_hsv(image_array)

    hsv_image += alteration
    hsv_image[hsv_image < 0] = 0
    hsv_image[hsv_image > 1] = 1

    modified_rgb_image = hsv_to_rgb(hsv_image)

    return modified_rgb_image


def rgb_image_to_rgb_0_1(image_array):
    one_value = image_array.flat[0]
    #print(one_value, type(one_value), isinstance(one_value, (int, numpy.integer)))
    if isinstance(one_value, (int, numpy.integer)):
        return image_array / 255
    else:
        return image_array
